#!/usr/bin/env python3
"""Validate and describe the minimal Wirehair2 benchmark contract.

This is intentionally a contract/ledger checker, not a campaign supervisor.
It keeps candidate experiments cheap while preventing an arm from improving
its denominator by omitting an inconvenient K, seed, or loss stratum.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple
import weakref


SCHEMA = "wirehair.wh2.benchmark-contract.v4"
FREEZE_SCHEMA = "wirehair.wh2.benchmark-freeze.v1"
LEGACY_RAW_FREEZE_SCHEMA = "wirehair.wh2.benchmark-freeze.v2"
RAW_FREEZE_SCHEMA = "wirehair.wh2.benchmark-freeze.v3"
REPAIR_MAP_SCHEMA = "wirehair.wh2.repair-map.v1"
TIMING_QUALIFICATION_MAP_SCHEMA = \
    "wirehair.wh2.timing-qualification-map.v1"
TIMING_QUALIFICATION_AUDIT_SCHEMA = \
    "wirehair.wh2.timing-qualification-audit.v1"
DEFAULT_CONTRACT = Path(__file__).with_name("wh2_benchmark_contract_v4.json")
MAX_JSON_LINE_BYTES = 64 * 1024
HEX64 = re.compile(r"0x[0-9a-f]{16}\Z")
HEX32 = re.compile(r"0x[0-9a-f]{8}\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
GIT_COMMIT = re.compile(r"[0-9a-f]{40}(?:[0-9a-f]{24})?\Z")
PROFILE_ID = 0x4b295bbb47f4f9c9
PROFILE_ENCODING_VERSION = 1
RAW_CONSTRUCTION_SEED_BASIS = "uniform-raw-v1"
PRODUCTION_CONSTRUCTION_SEED_BASIS = "production-profile-v1"
NOT_APPLICABLE_CONSTRUCTION_SEED_BASIS = "not-applicable"
RAW_SEED_SCHEDULE_SHA256 = (
    "90a98a3db207852dabdf5fb27573ef48bce52e0228cee4e291d96fa44ed509a7"
)
RAW_PRECODE_BASE_SEED = 0x487468302aad7105
RAW_PRECODE_ATTEMPT_STRIDE = 0x9e3779b97f4a7c15
RAW_PACKET_BASE_SEED = 0x4ec72102
RAW_PACKET_ATTEMPT_STRIDE = 0x9e3779b9
LEGACY_RAW_ARM_DESCRIPTOR_SHA256S = frozenset((
    # Legacy v3 raw descriptors included the seed schedule inside the
    # descriptor. They remain denylisted so removing the raw arm-name prefix
    # cannot downgrade an old receipt into generic v1 evidence.
    "0550e0ed0c62d5491ff6915652fd96ed25f3c7782462da8c551636ec2e0294dd",
    "80f57b83eac9b8e3a19d8235cc067e39990980510e46588ddefa16f9561e1c38",
    "0eb3aef0602b5e7de15c822de84a5dbfc5dfdd99b76fbfd41538f7a13248c3a5",
    "2dc244661b3b073569319377ee3e55333a82ddad7bd328e1b0fef67395174614",
))
RAW_STRUCTURE_BY_DESCRIPTOR = {
    # Raw v2 keeps structure identity separate from its seed schedule.  Each
    # closed descriptor hashes its canonical arm name and equation transform,
    # so the validator must bind both that name and the realized fixed
    # geometry before accepting a recomputed construction receipt.
    "91d7c1a558e1cf93b002fcf2062b7657d301faca03972215495bdf2429499e90": {
        "arm": "wirehair2_raw_d12_h11_periodic",
        "binary_dense_rows": 12, "gf256_heavy_rows": 11,
        "heavy_family": "periodic-cauchy", "mix_count": 3,
        "dense_anchor_layout": "disabled",
    },
    "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11": {
        "arm": "wirehair2_raw_d12_h12_periodic",
        "binary_dense_rows": 12, "gf256_heavy_rows": 12,
        "heavy_family": "periodic-cauchy", "mix_count": 3,
        "dense_anchor_layout": "disabled",
    },
    "7c7889747a97ac160726b807fb03349344d49d4bec84c9e8220aa4689b00d2ca": {
        "arm": "wirehair2_raw_d12_h13_periodic",
        "binary_dense_rows": 12, "gf256_heavy_rows": 13,
        "heavy_family": "periodic-cauchy", "mix_count": 3,
        "dense_anchor_layout": "disabled",
    },
    "c70e0f57bb8d7783fa29b0decbed5da5058a8eb532d57d540f72108e114f091a": {
        "arm": "wirehair2_raw_d13_h12_periodic",
        "binary_dense_rows": 13, "gf256_heavy_rows": 12,
        "heavy_family": "periodic-cauchy", "mix_count": 3,
        "dense_anchor_layout": "disabled",
    },
    "9527f200ad38c7eec6502b2f768fdd67b92787fb227eed3d7616274ffc2df388": {
        "arm": "wirehair2_dense_two07_basis_v1",
        "binary_dense_rows": 12, "gf256_heavy_rows": 12,
        "heavy_family": "periodic-cauchy", "mix_count": 3,
        "dense_anchor_layout": "two07",
    },
}
RAW_ARM_DESCRIPTOR_SHA256S = LEGACY_RAW_ARM_DESCRIPTOR_SHA256S | frozenset(
    RAW_STRUCTURE_BY_DESCRIPTOR)
MEASURED_BATCH_BLOCK_TARGET = 65536
INVOCATIONS_PER_SLOT_FORMULA = (
    "max(2,ceil(measured_batch_block_target/K))"
)
TIMING_INTERLEAVE_POLICY = "self-counterbalanced-repeat-major-v1"
TIMING_LOSS_RETRY_STRIDE = 0x9e3779b97f4a7c15
TIMING_LOSS_RETRY_ATTEMPTS = 256
TIMING_QUALIFICATION_ENTRY_KIND = \
    "loss_retry_offset_indexed_by_base_timing_cell_ordinal"
SCORED_OUTCOMES = ("success", "need_more_at_cap", "construct_failed", "unsupported")
EXPECTED_OUTCOMES = {
    "success": "decoded_extra is an integer in [0,overhead_cap]",
    "need_more_at_cap": "explicit failure charged overhead_cap",
    "construct_failed": "explicit failure charged overhead_cap",
    "unsupported": "explicit failure charged overhead_cap",
    "fatal_or_internal": "abort the arm; never score",
}
EXPECTED_BANDS = (
    ("2-100", 2, 100), ("101-1000", 101, 1000),
    ("1001-5000", 1001, 5000), ("5001-10000", 5001, 10000),
    ("10001-20000", 10001, 20000), ("20001-64000", 20001, 64000),
)
EXPECTED_SHORT_K = (
    2, 3, 4, 5, 6, 8, 16, 32, 64, 100,
    101, 128, 256, 512, 513, 1000,
    1001, 2048, 4096, 5000,
    5001, 8192, 10000,
    10001, 16384, 20000,
    20001, 32768, 49152, 64000,
)
_RAW_TINY_STAIRCASE = (
    0, 0, 2, 3, 3, 5, 6, 6, 6, 7, 9, 10, 10, 10, 12, 14,
    13, 14, 12, 12, 15, 16, 21, 14, 14, 13, 18, 21, 22, 21, 13, 22,
    13, 24, 14, 17, 16, 24, 30, 26, 24, 18, 15, 15, 24, 18, 21, 17,
    14, 16, 21, 18, 17, 22, 25, 20, 17, 18, 21, 18, 23, 20, 19, 23,
    19,
)
_RAW_STAIRCASE_POINTS = (
    (2048, 52), (2618, 54), (2826, 60), (3725, 62),
    (3962, 67), (4277, 65), (4547, 60), (5065, 64),
    (5224, 76), (5642, 66), (5909, 71), (6285, 76),
    (6583, 66), (6895, 72), (7448, 69), (7682, 76),
    (8046, 78), (8558, 76), (8963, 73), (9389, 81),
    (10143, 86), (11129, 94), (12593, 99), (12988, 105),
    (14032, 108), (14473, 114), (15397, 110), (16636, 113),
    (17698, 118), (18828, 123), (19420, 127), (20343, 136),
    (21979, 139), (23024, 150), (24119, 156), (25659, 162),
    (27298, 173), (29042, 176), (30898, 183), (31870, 190),
    (33906, 200), (35519, 211), (37208, 220), (38978, 234),
    (40205, 253), (42776, 297), (44122, 320), (45511, 336),
    (46944, 357), (48421, 373), (49177, 376), (50725, 380),
    (52321, 391), (53968, 388), (54811, 382), (54811, 382),
    (55667, 372), (57419, 362), (58316, 356), (60152, 347),
    (61091, 337), (62045, 334), (63014, 340), (64000, 345),
)


def _raw_staircase_for_K(K: int) -> int:
    """Exact Python port of GetDenseCount(), used by raw WH2 structure."""
    if type(K) is not int or not 2 <= K <= 64000:
        fail("raw staircase K is outside [2,64000]")
    if K < len(_RAW_TINY_STAIRCASE):
        return _RAW_TINY_STAIRCASE[K]
    if K < 2048:
        if K <= 500:
            low_K, high_K, low_count, high_count = 64, 500, 26, 35
        elif K <= 1000:
            low_K, high_K, low_count, high_count = 500, 1000, 35, 48
        else:
            low_K, high_K, low_count, high_count = 1000, 2048, 48, 62
    else:
        low = 0
        high = len(_RAW_STAIRCASE_POINTS) - 1
        while True:
            middle = (high + low) // 2
            if middle == low:
                break
            if K > _RAW_STAIRCASE_POINTS[middle][0]:
                low = middle
            else:
                high = middle
        low_K, low_count = _RAW_STAIRCASE_POINTS[low]
        high_K, high_count = _RAW_STAIRCASE_POINTS[low + 1]
    numerator = (K - low_K) * (high_count - low_count)
    denominator = high_K - low_K
    quotient = numerator // denominator if numerator >= 0 else \
        -((-numerator) // denominator)
    count = low_count + quotient
    remainder = count % 4
    if remainder == 0:
        count += 2
    elif remainder == 1:
        count += 1
    elif remainder == 3:
        count += 3
    return count


RAW_STAIRCASE_BY_K = {
    K: _raw_staircase_for_K(K) for K in EXPECTED_SHORT_K
}
EXPECTED_TIMING_SHORT_K = (
    8, 32, 100, 128, 512, 1000,
    2048, 5000, 8192, 20000, 32768, 64000,
)
EXPECTED_RAW_BASE_ATTEMPTS = (0, 1, 2)
EXPECTED_TRAINING_LOSS_ROOTS = (
    "0xd1b54a32d192ed03", "0x94d049bb133111eb", "0x8538ecb5bd456ea3",
)
EXPECTED_VALIDATION_LOSS_ROOTS = (
    "0xc0ac29b7c97c50dd", "0x3f84d5b5b5470917", "0x9216d5d98979fb1b",
)
EXPECTED_STRATA_SETS = {
    "development": (
        ("iid", 100000),
        ("burst", 500000),
        ("adversarial", 500000),
        ("repair-only", 500000),
    ),
    "hard": (
        ("burst", 500000),
        ("adversarial", 500000),
        ("repair-only", 500000),
    ),
}
EXPECTED_RECOVERY_DOMAINS = {
    "development": ("short", (2,), "raw_paired_training", "development", 360),
    "final_raw": ("all", (2,), "raw_paired_training", "hard", 575991),
    "final_repaired": ("all", (2,), "production_training", "hard", 575991),
    "final_validation": ("all", (2,), "production_validation", "hard", 575991),
    "cross_width_validation": (
        "short", (64, 256, 1280, 4096), "production_validation",
        "development", 1440,
    ),
}
EXPECTED_PACKET_TRACE = {
    "loss_rate_conversion": "double(loss_ppm) / 1000000.0",
    "trace_seed": (
        "loss_seed xor (K * 0x9e3779b97f4a7c15) xor "
        "(block_bytes * 0xbf58476d1ce4e5b9) mod 2^64; overhead salt and "
        "local trial salt are zero"
    ),
    "rng": "SplitMix64; Unit=(Next() >> 11) / 2^53; drop iff Unit < loss_rate",
    "iid": "candidate packet IDs 0,1,...; use the unsalted trace-seed RNG",
    "burst": (
        "candidate packet IDs 0,1,...; use RNG seed trace_seed xor 0x10fade; "
        "eight-drop bursts start on an idle candidate with probability "
        "loss_rate/(8-7*loss_rate)"
    ),
    "adversarial": (
        "candidate packet IDs UINT32_MAX-2*i modulo 2^32; use RNG seed "
        "trace_seed xor 0x10fade and IID drops"
    ),
    "repair-only": (
        "candidate packet IDs K+i modulo 2^32; use RNG seed trace_seed xor "
        "0x10fade and IID drops"
    ),
    "nested_overhead": (
        "generate one K+4 delivered-ID prefix; threshold h replays exactly "
        "its first K+h IDs"
    ),
    "candidate_limit": "256*(K+4)+65536 attempted candidates",
    "trace_sha256": (
        "SHA-256 of the K+4 delivered IDs encoded consecutively as unsigned "
        "32-bit little-endian words"
    ),
}
EXPECTED_PHASE_SEED_POLICY = {
    "development": "one_base_attempt_no_retry",
    "final_raw": "one_base_attempt_no_retry",
    "final_repaired": "frozen_lowest_training_success_attempt",
    "final_validation": "replay_frozen_training_repair_map",
    "cross_width_validation": "replay_frozen_training_repair_map",
}
EXPECTED_REPAIR_RULE = (
    "for each K, choose the lowest production retry offset in [0,255] that "
    "decodes at overhead 0 on all three hard schedules crossed with all three "
    "training loss roots; freeze the complete map and SHA-256 before "
    "validation; no manual exceptions"
)
EXPECTED_ATTEMPT_DERIVATION = (
    "retry offset r selects public seed_attempt=(base_seed_attempt+r) modulo "
    "256; that public "
    "attempt adds attempt*0x9e3779b97f4a7c15 modulo 2^64 to the certified "
    "GF(256) precode seed and attempt*0x9e3779b9 modulo 2^32 to the packet "
    "peel seed; raw phases use r=0 exactly"
)
EXPECTED_TIMING_QUALIFICATION = {
    "map_schema": TIMING_QUALIFICATION_MAP_SCHEMA,
    "audit_schema": TIMING_QUALIFICATION_AUDIT_SCHEMA,
    "entry_kind": TIMING_QUALIFICATION_ENTRY_KIND,
    "max_attempts": TIMING_LOSS_RETRY_ATTEMPTS,
    "retry_seed_derivation": (
        "loss retry offset r in [0,255] selects loss_seed="
        "(base_loss_seed+r*0x9e3779b97f4a7c15) modulo 2^64; r=0 "
        "preserves the base trace"
    ),
    "selection_rule": (
        "for each base timing cell choose the lowest retry offset where "
        "wirehair2_head byte-recovers from the exact first K+4 delivered "
        "IDs and wirehair1 byte-recovers by K+256; no cell may be omitted"
    ),
    "retryable_outcome": (
        "only need_more_at_bound is retryable; construction, unsupported, "
        "fatal, internal, allocation, terminal, and byte-mismatch outcomes "
        "abort qualification"
    ),
    "candidate_policy": (
        "qualification runs exactly the two mandatory controls before any "
        "candidate timing; candidate arms, outcomes, and timings are "
        "forbidden and can never select or retry a loss trace"
    ),
    "controls": [
        {
            "arm": "wirehair2_head",
            "scope": "decoder_solve",
            "success": "byte-exact recovery at decoded_extra=4",
        },
        {
            "arm": "wirehair1",
            "scope": "receive_to_success",
            "success": "byte-exact recovery at decoded_extra in [0,256]",
        },
    ],
}
EXPECTED_TIMING_POLICY = {
    "development_wall_time_seconds_per_candidate": 7200,
    "loss_ppm": 100000,
    "schedule": "iid",
    "fixed_received_overhead": 4,
    "receive_overhead_cap": 256,
    "seed_derivation": (
        "construction seed is selected by seed_mode and unchanged; "
        "base_loss_seed=SplitMix64(loss root xor (replicate * "
        "0x9e3779b97f4a7c15 mod 2^64)); realized loss_seed is selected "
        "only by the frozen timing qualification map"
    ),
    "source_derivation": (
        "source bytes are derived only from canonical base_cell_key JSON "
        "under wirehair.wh2.source.v1; base_cell_sha256 authenticates that "
        "JSON and loss retry fields never enter the source identity"
    ),
    "confidence": (
        "two-sided Student-t 95% interval over independent-round mean paired "
        "log seconds"
    ),
    "practical_margin_ppm": 20000,
    "effective_floor": "log1p(0.02) after both A/A intervals pass",
    "aa_repeatability_rule": (
        "each A/A 95% interval must lie strictly inside "
        "[-log1p(0.02),+log1p(0.02)]; otherwise the comparison is "
        "non-selectable"
    ),
    "order": "self-counterbalanced primary ABBA/BAAB plus opposite order",
    "cache_state": "warm",
    "production_policy": (
        "final timing replays the frozen production repair map; each "
        "Wirehair2 timing receipt binds its construction map SHA-256 and "
        "realized attempt; every timing arm binds the phase qualification map"
    ),
    "eligibility": (
        "the frozen candidate-blind qualification map guarantees mandatory "
        "control success; any candidate failure invalidates its timing panel "
        "without changing the trace, retrying, dropping a row, or forming an "
        "observed common-success intersection"
    ),
    "scope_protocol": (
        "one qualified arm-free K+256 delivered-ID trace is frozen per cell; "
        "decoder_solve uses a fresh decoder and its identical first K+4 "
        "received prefix; encoder_init_plus_first_K_symbols encodes exactly "
        "IDs 0..K-1; receive_to_success replays the trace and includes all "
        "feed and final recovery work through each arm's first success, with "
        "need_more_at_cap after K+256"
    ),
}
EXPECTED_PANEL_PROTOCOL = {
    "control_aa": [
        {"arm": "wirehair2_head", "scope": "decoder_solve"},
        {"arm": "wirehair1", "scope": "receive_to_success"},
        {"arm": "wirehair2_head",
         "scope": "encoder_init_plus_first_K_symbols"},
        {"arm": "wirehair1",
         "scope": "encoder_init_plus_first_K_symbols"},
    ],
    "candidate_aa_scopes": [
        "decoder_solve", "receive_to_success",
        "encoder_init_plus_first_K_symbols",
    ],
    "candidate_ab": [
        {"control": "wirehair2_head", "scope": "decoder_solve"},
        {"control": "wirehair1", "scope": "receive_to_success"},
        {"control": "wirehair2_head",
         "scope": "encoder_init_plus_first_K_symbols"},
        {"control": "wirehair1",
         "scope": "encoder_init_plus_first_K_symbols"},
    ],
    "slots_per_panel": 8,
    "warmups_per_logical_side": 1,
    "measured_batch_block_target": MEASURED_BATCH_BLOCK_TARGET,
    "invocations_per_slot": INVOCATIONS_PER_SLOT_FORMULA,
    "interleave_policy": TIMING_INTERLEAVE_POLICY,
    "order_assignment": (
        "phase_bit=low bit of SHA256(canonical panel key); primary order is "
        "ABBA iff lane parity equals phase_bit, otherwise BAAB; slots 0..3 "
        "use primary order and slots 4..7 use its opposite"
    ),
    "batch_split": (
        "n=invocations_per_slot=max(2,ceil(measured_batch_block_target/K)); "
        "each primary-order slot executes ceil(n/2) calls and each "
        "opposite-order slot executes floor(n/2) calls"
    ),
    "within_block_traversal": "repeat-major",
}
EXPECTED_TIMING_EXECUTION_GEOMETRY = {
    "timing_worker_count": 8,
    "jobs_per_wave": 8,
    "cohort_identity": (
        "all timing cell fields except replicate, base_loss_seed, "
        "base_cell_sha256, loss_retry_offset, and loss_seed, plus panel"
    ),
    "cohort_order": (
        "block_bytes order, then K-set order, then frozen panel ordinal"
    ),
    "round_partition": "independent_round=floor(replicate/8); lane=replicate%8",
    "round_traversal": (
        "independent round outermost; traverse the complete frozen cohort "
        "domain in cohort order; execute one eight-lane wave per cohort"
    ),
    "worker_slot": "(lane+cohort_index+independent_round)%8",
    "barrier": "after_each_eight-lane_cohort_wave",
    "cohort_count_formula": (
        "(expected_cells/paired_repetitions)*"
        "len(timing_panels(contract,frozen_arms))"
    ),
    "topology_policy": (
        "eight singleton-pinned workers on distinct physical cores; maximize "
        "coverage of available L3/LLC groups, then fill deterministically; "
        "distinct L3/LLC groups are preferred but not required"
    ),
}
EXPECTED_TIMING_STATISTICS = {
    "lane_log_ratio": (
        "average the left/right-oriented four-slot log contrast from primary "
        "slots 0..3 and the opposite-order contrast from slots 4..7"
    ),
    "independent_unit": (
        "one independent-round mean formed from all eight lane replicate "
        "means"
    ),
    "aggregate_weighting": (
        "within each lane, equal weight for every frozen K in the reported "
        "band and width; aggregate timing equal-weights every frozen K x "
        "width cell; then equal-weight all eight lanes in the round"
    ),
    "variance": "unbiased sample variance over independent-round means",
    "t_critical_by_independent_rounds": {
        "12": 2.200985160082949,
    },
    "block_phase_difference": (
        "primary four-slot log contrast minus opposite-order four-slot log "
        "contrast; this diagnoses first-block versus second-block drift and "
        "is summarized over the same 12 independent round means"
    ),
    "faster": "A/B upper 95% bound < -effective_floor",
    "noninferior": "A/B upper 95% bound < effective_floor",
    "resolved_slower": "A/B lower 95% bound > effective_floor",
    "strict_inequalities": True,
}
EXPECTED_INVALID_TIMER_POLICY = (
    "let n=invocations_per_slot; each successful primary-order slot sums "
    "exactly ceil(n/2) fresh exact-scope invocation-reported durations and "
    "each successful opposite-order slot sums exactly floor(n/2), producing "
    "eight integer elapsed_ns values in [1,2^63-1]; any non-success requires "
    "eight nulls and makes the complete comparison panel non-selectable; rows "
    "are never dropped"
)
EXPECTED_SELECTION_POLICY = {
    "controls": ["wirehair2_head", "wirehair1"],
    "development_architecture_roles": {
        "descriptive_controls": ["wirehair2_head", "wirehair1"],
        "recovery_reference": "wirehair2_raw_d12_h12_periodic",
        "recovery_candidates": ["wirehair2_dense_two07_basis_v1"],
        "timing_proxy": "wirehair2_head",
        "timing_candidates": ["wirehair2_dense_two07_basis_v1"],
    },
    "development_recovery_eligibility": (
        "candidate eligibility is gated only against the uniform-raw D12/H12 "
        "recovery reference; production WH2 and Wirehair1 are descriptive "
        "non-vetoing controls"
    ),
    "development_timing_eligibility": (
        "time only the frozen Two07 candidate against the production-head "
        "D12 structure proxy; the raw D12 recovery reference has no timing "
        "panel"
    ),
    "development_timing_proxy_policy": (
        "the witness proves only D12-disabled non-seed structure equivalence "
        "under development attempt-0 production timing seeds, separately "
        "binds the distinct uniform-raw recovery seed schedule, and requires "
        "a new witness for any other attempt semantics"
    ),
    "noninferiority_margin_ppm": 1000,
    "architecture_failure_equivalence_ppm": 1000,
    "raw_weak_seed_definition": (
        "a (K,base_seed_attempt,block_bytes) with any overhead-0 failure"
    ),
    "raw_individual_weak_seeds_are_vetoes": False,
    "raw_repairs_and_introductions_are_descriptive": True,
    "architecture_order": [
        "eligible aggregate/band/stratum overhead-tail noninferiority",
        "form the recovery-equivalent set within 0.1 percentage point of the "
        "eligible minimum aggregate overhead-0 failure count",
        "fastest decoder_solve within that recovery-equivalent set",
        "lexicographic failures at overhead 1,2,4",
        "stable arm identifier",
    ],
    "development_resolved_slowdown_rejects": True,
    "final_decoder_solve_vs_wirehair2": (
        "candidate/wirehair2_head paired-log 95% upper bound is strictly "
        "below negative effective floor in every band x payload width"
    ),
    "final_receive_to_success_vs_wirehair1": (
        "candidate/wirehair1 paired-log 95% upper bound is strictly below "
        "negative effective floor in every band x payload width"
    ),
    "final_encoder_rule": (
        "candidate/control paired-log 95% upper bound is strictly below "
        "effective floor for both controls in every band x payload width and "
        "strictly below negative effective floor for both controls under "
        "equal-cell aggregate weighting"
    ),
    "final_raw_failure_rule": (
        "100 * failures <= cells overall and in every band x loss/schedule "
        "stratum at overhead 0"
    ),
    "final_repaired_rule": (
        "zero overhead-0 weak (K,production_seed) units on the training census"
    ),
    "final_validation_rule": (
        "zero overhead-0 weak K units on disjoint loss roots with zero "
        "unsupported, construction, fatal, or internal errors; the 1% rule "
        "is therefore implied"
    ),
}
EXPECTED_EVIDENCE_POLICY = {
    "freeze_before_results": [
        "all freezes: source_git_commit", "all freezes: contract_sha256",
        "all freezes: arm_roster_sha256", "all freezes: domain_sha256",
        "all freezes: trace_manifest_sha256",
        "all freezes per arm: binary_sha256",
        "all freezes per arm: arm_descriptor_sha256",
        "all freezes per arm: repair_map_sha256",
        "all freezes: repair_training_trace_manifest_sha256",
        "all freezes: commands", "all freezes: cpu_affinity",
        "all freezes: host_identity",
        "development raw-v3 per arm: construction_seed_basis",
        "development raw-v3 per arm: seed_schedule_sha256",
        "development raw-v3 per arm: dense_anchor_layout",
        "development raw-v3: architecture_roles",
        "development raw-v3: timing_proxy_witness_sha256",
        "development raw-v3: work_rank_summary_sha256",
        "development raw-v3: work_rank_result_stream_sha256",
        "development raw-v3: work_rank_domain_sha256",
    ],
    "freeze_manifest_schema": (
        "phase-scoped canonical JSON: benchmark-freeze.v3 is mandatory for "
        "development raw recovery; benchmark-freeze.v1 remains mandatory "
        "for timing and non-development production phases; read-only legacy "
        "raw-v2 development evidence cannot authenticate dense-anchor arms; "
        "each accepted manifest is the sole arm-roster and artifact authority"
    ),
    "trace_manifest_schema": (
        "recovery uses wirehair.wh2.trace-manifest.v1; timing uses v2 with "
        "a hash header binding the base domain, qualified domain, and "
        "qualification-map SHA-256; both use ordered canonical JSONL rows of "
        "ordinal,cell_sha256,trace_sha256 and permit duplicate trace hashes"
    ),
    "repair_map_schema": (
        "wirehair.wh2.repair-map.v1 canonical JSON object with 63999 retry "
        "offsets indexed by K-2 and bound to training traces"
    ),
    "timing_qualification_schema": (
        "one phase-scoped wirehair.wh2.timing-qualification-map.v1 canonical "
        "JSON object binds an ordered candidate-blind audit and one lowest "
        "uint8 loss retry offset per base timing-cell ordinal"
    ),
    "timing_qualification_binding": (
        "the timing freeze directly binds the qualified timing domain and "
        "timing trace-manifest v2 SHA-256; that v2 hash header transitively "
        "binds the base domain, qualified domain, and validated qualification-"
        "map SHA-256, whose canonical object binds the complete audit SHA-256"
    ),
    "canonical_hashing": (
        "parse with duplicate-key rejection; accept LF or CRLF input; hash "
        "sorted compact logical JSON with LF and explicit schema domain separation"
    ),
    "common_cell_policy": (
        "every arm has exactly one result for every qualified arm-free key; "
        "source identity is derived only from its authenticated base cell"
    ),
    "invalid_domain_mutations": [
        "missing", "duplicate", "extra", "seed_swapped", "partial",
        "wrong_band", "wrong_loss", "incomplete_roster", "domain_hash_drift",
        "trace_drift", "repair_map_drift", "construction_attempt_drift",
        "qualification_map_drift", "qualification_audit_drift",
        "loss_retry_drift", "source_identity_drift",
        "timing_panel_drift", "timing_order_drift", "timing_batch_drift",
        "timing_geometry_drift", "raw_schema_downgrade",
        "dense_anchor_layout_drift", "architecture_role_drift",
        "timing_proxy_witness_drift", "work_rank_binding_drift",
    ],
    "selectable_excluded_cells": 0,
    "unsupported_policy": (
        "only an explicitly measured routed composite is selectable"
    ),
    "thermal_policy": (
        "run under the existing sole CPU/DIMM/EDAC sampler; abort on its "
        "configured thermal or EDAC violation"
    ),
    "final_continuity_policy": (
        "one joint continuity check must bind final_raw, final_repaired, "
        "final_validation, cross_width_validation, and final timing to the "
        "same source, host, roster, codec kinds, binaries, descriptors, "
        "training trace, and production repair maps, and to the exact "
        "development architecture-selection receipt"
    ),
}
LEDGER_FIELDS = frozenset((
    "arm", "phase", "band", "K", "block_bytes", "loss_ppm", "schedule",
    "trial", "base_seed_attempt", "loss_seed", "overhead_cap", "outcome",
    "decoded_extra", "cell_sha256", "trace_sha256", "binary_sha256",
    "arm_descriptor_sha256", "construction_attempt",
    "realized_construction_sha256", "repair_map_sha256",
))
LEGACY_RAW_CONSTRUCTION_FIELDS = (
    "construction_seed_basis", "seed_schedule_sha256",
    "precode_attempt", "packet_attempt", "effective_precode_seed",
    "effective_packet_seed", "staircase", "binary_dense_rows",
    "gf256_heavy_rows", "source_hits", "dense_identity_corner",
    "heavy_family", "mix_count",
)
RAW_CONSTRUCTION_FIELDS = LEGACY_RAW_CONSTRUCTION_FIELDS + (
    "dense_anchor_layout",
)
LEGACY_RAW_RECOVERY_RECORD_FIELDS = LEDGER_FIELDS | frozenset(
    LEGACY_RAW_CONSTRUCTION_FIELDS)
RAW_RECOVERY_RECORD_FIELDS = LEDGER_FIELDS | frozenset(
    RAW_CONSTRUCTION_FIELDS)
TIMING_RECEIPT_FIELDS = frozenset((
    "phase", "band", "K", "block_bytes", "loss_ppm", "schedule",
    "replicate", "base_seed_attempt", "base_loss_seed", "base_cell_sha256",
    "loss_retry_offset", "loss_seed", "fixed_received_overhead", "receive_overhead_cap",
    "invocations_per_slot", "interleave_policy",
    "panel_kind", "scope", "left_arm", "right_arm", "order",
    "left_outcome", "right_outcome", "left_decoded_extra",
    "right_decoded_extra", "elapsed_ns", "cell_sha256", "trace_sha256",
    "left_binary_sha256", "right_binary_sha256",
    "left_arm_descriptor_sha256", "right_arm_descriptor_sha256",
    "left_construction_attempt", "right_construction_attempt",
    "left_realized_construction_sha256",
    "right_realized_construction_sha256", "left_repair_map_sha256",
    "right_repair_map_sha256",
))
TRACE_FIELDS = frozenset((
    "ordinal", "cell_sha256", "trace_sha256",
))
FREEZE_FIELDS = frozenset((
    "schema", "contract_sha256", "evidence_kind", "phase",
    "domain_sha256", "source_git_commit", "arm_roster",
    "arm_roster_sha256", "trace_manifest_sha256",
    "repair_training_trace_manifest_sha256", "commands",
    "cpu_affinity", "host_identity", "arms",
))
FREEZE_ARM_FIELDS = frozenset((
    "arm", "codec", "binary_sha256", "arm_descriptor_sha256",
    "construction_policy", "repair_map_sha256",
))
RAW_FREEZE_ARM_FIELDS = FREEZE_ARM_FIELDS | frozenset((
    "construction_seed_basis", "seed_schedule_sha256",
))
RAW_V3_FREEZE_ARM_FIELDS = RAW_FREEZE_ARM_FIELDS | frozenset((
    "dense_anchor_layout",
))
RAW_V3_FREEZE_FIELDS = FREEZE_FIELDS | frozenset((
    "architecture_roles", "timing_proxy_witness_sha256",
    "work_rank_summary_sha256", "work_rank_result_stream_sha256",
    "work_rank_domain_sha256",
))
RAW_ARCHITECTURE_ROLE_FIELDS = frozenset((
    "descriptive_controls", "recovery_reference", "recovery_candidates",
    "timing_proxy", "timing_candidates",
))
EXPECTED_RAW_ARCHITECTURE_ROLES = {
    "descriptive_controls": ["wirehair2_head", "wirehair1"],
    "recovery_reference": "wirehair2_raw_d12_h12_periodic",
    "recovery_candidates": ["wirehair2_dense_two07_basis_v1"],
    "timing_proxy": "wirehair2_head",
    "timing_candidates": ["wirehair2_dense_two07_basis_v1"],
}
REPAIR_MAP_FIELDS = frozenset((
    "schema", "contract_sha256", "training_domain_sha256", "arm",
    "source_git_commit", "binary_sha256", "arm_descriptor_sha256",
    "production_base_seed_attempt", "entry_kind", "attempt_derivation",
    "repair_rule", "training_trace_manifest_sha256", "retry_offsets",
))
TIMING_QUALIFICATION_MAP_FIELDS = frozenset((
    "schema", "contract_sha256", "phase", "source_git_commit",
    "base_domain_sha256", "qualified_domain_sha256", "entry_kind",
    "controls", "qualification_audit_sha256",
    "selected_trace_roster_sha256", "retry_offsets",
))
TIMING_QUALIFICATION_CONTROL_FIELDS = frozenset((
    "arm", "scope", "binary_sha256", "arm_descriptor_sha256",
    "construction_policy", "repair_map_sha256",
))
TIMING_QUALIFICATION_CONTROL_ORDER = (
    "arm", "scope", "binary_sha256", "arm_descriptor_sha256",
    "construction_policy", "repair_map_sha256",
)
TIMING_QUALIFICATION_AUDIT_FIELDS = frozenset((
    "ordinal", "base_cell_sha256", "loss_retry_offset", "loss_seed",
    "trace_sha256", "wirehair2_head_outcome",
    "wirehair2_head_decoded_extra", "wirehair1_outcome",
    "wirehair1_decoded_extra",
))
FINAL_FREEZE_KEYS = frozenset((
    ("recovery", "final_raw"),
    ("recovery", "final_repaired"),
    ("recovery", "final_validation"),
    ("recovery", "cross_width_validation"),
    ("timing", "final"),
))
SELECTION_FIELDS = frozenset((
    "schema", "contract_sha256", "recovery_domain_sha256",
    "timing_base_domain_sha256", "timing_domain_sha256",
    "timing_qualification_map_sha256", "recovery_freeze_manifest_sha256",
    "timing_freeze_manifest_sha256", "architecture_artifact_sha256",
    "recovery_cells_per_arm", "timing_rows", "candidate_roster",
    "eligible_candidates", "eligible_overhead0_failures",
    "minimum_overhead0_failures", "recovery_equivalence_allowance",
    "recovery_equivalent_candidates", "ranking", "selected_arm",
    "selected_codec", "selected_arm_descriptor_sha256",
    "selected_architecture_sha256", "selection_sha256",
))


class ContractError(RuntimeError):
    """The frozen benchmark contract or a result ledger is invalid."""


TIMING_QUALIFICATION_FIELDS = (
    "contract_sha256", "phase", "map_sha256", "base_domain_sha256",
    "qualified_domain_sha256", "source_git_commit", "controls",
    "qualification_audit_sha256", "retry_offsets",
    "selected_trace_roster_sha256", "selected_trace_sha256s",
)


class TimingQualification:
    """Opaque handle for a strictly loaded phase timing qualification."""

    __slots__ = TIMING_QUALIFICATION_FIELDS + ("__weakref__",)

    def __new__(cls, *_args: Any, **_kwargs: Any) -> "TimingQualification":
        fail("TimingQualification objects must come from the strict loader")
        raise AssertionError("unreachable")

    def __setattr__(self, _name: str, _value: Any) -> None:
        raise AttributeError("TimingQualification is immutable")


def _timing_qualification_registry() -> Tuple[Any, Any]:
    """Return a loader-only sealer and an out-of-band provenance verifier."""
    snapshots: Any = weakref.WeakKeyDictionary()
    fields = TIMING_QUALIFICATION_FIELDS

    def snapshot(value: TimingQualification) -> Tuple[Any, ...]:
        return tuple(getattr(value, field)
                     for field in fields)

    def seal(values: Mapping[str, Any]) -> TimingQualification:
        if set(values) != set(fields):
            fail("internal timing qualification seal has the wrong schema")
        value = object.__new__(TimingQualification)
        for field in fields:
            object.__setattr__(value, field, values[field])
        snapshots[value] = snapshot(value)
        return value

    def verify(value: Any) -> bool:
        if type(value) is not TimingQualification:
            return False
        try:
            expected = snapshots.get(value)
            return expected is not None and expected == snapshot(value)
        except (AttributeError, TypeError):
            return False

    return seal, verify


_seal_timing_qualification, _timing_qualification_is_registered = \
    _timing_qualification_registry()


def fail(message: str) -> None:
    raise ContractError(message)


def _object_no_duplicates(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            fail("duplicate JSON key: " + key)
        result[key] = value
    return result


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
    )


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def repair_map_sha256(value: Mapping[str, Any]) -> str:
    digest = hashlib.sha256()
    digest.update(b"wirehair.wh2.repair-map.v1\0")
    digest.update(canonical_json(value).encode("utf-8"))
    return digest.hexdigest()


def timing_qualification_map_sha256(value: Mapping[str, Any]) -> str:
    digest = hashlib.sha256()
    digest.update(TIMING_QUALIFICATION_MAP_SCHEMA.encode("ascii") + b"\0")
    digest.update(canonical_json(value).encode("utf-8"))
    return digest.hexdigest()


def timing_selected_trace_roster_sha256(values: Sequence[str]) -> str:
    if (not isinstance(values, Sequence) or
            isinstance(values, (str, bytes)) or
            any(not isinstance(value, str) or SHA256.fullmatch(value) is None
                for value in values)):
        fail("selected timing traces must be lowercase SHA-256 values")
    digest = hashlib.sha256()
    digest.update(b"wirehair.wh2.timing-selected-trace-roster.v1\0")
    digest.update(canonical_json(list(values)).encode("utf-8"))
    return digest.hexdigest()


def _timing_qualification_audit_hasher(
        contract: Mapping[str, Any], phase: str) -> Any:
    domains = contract.get("timing", {}).get("domains", {}) \
        if isinstance(contract, Mapping) else {}
    domain = domains.get(phase) if isinstance(domains, Mapping) else None
    if not isinstance(domain, Mapping):
        fail("unknown timing qualification phase: " + str(phase))
    base_domain_sha256 = domain.get("base_domain_sha256")
    if (not isinstance(base_domain_sha256, str) or
            SHA256.fullmatch(base_domain_sha256) is None):
        fail("timing qualification phase lacks a base-domain SHA-256")
    digest = hashlib.sha256()
    digest.update(TIMING_QUALIFICATION_AUDIT_SCHEMA.encode("ascii") + b"\0")
    digest.update(bytes.fromhex(contract_sha256(contract)))
    digest.update(bytes.fromhex(base_domain_sha256))
    digest.update(phase.encode("utf-8"))
    digest.update(b"\0")
    return digest


def timing_qualification_audit_sha256(
        contract: Mapping[str, Any], phase: str, path: Path) -> str:
    digest = _timing_qualification_audit_hasher(contract, phase)
    for value in _parse_canonical_jsonl(path, "timing qualification audit"):
        digest.update(canonical_json(value).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def arm_roster_sha256(arms: Sequence[str]) -> str:
    digest = hashlib.sha256()
    digest.update(b"wirehair.wh2.arm-roster.v1\0")
    digest.update(canonical_json(list(arms)).encode("utf-8"))
    return digest.hexdigest()


def freeze_manifest_sha256(value: Mapping[str, Any]) -> str:
    manifest = {
        key: item for key, item in value.items() if key != "arms_by_name"
    }
    schema = manifest.get("schema")
    if schema not in (
            FREEZE_SCHEMA, LEGACY_RAW_FREEZE_SCHEMA, RAW_FREEZE_SCHEMA):
        fail("cannot hash an unknown freeze-manifest schema")
    digest = hashlib.sha256()
    digest.update(schema.encode("ascii") + b"\0")
    digest.update(canonical_json(manifest).encode("utf-8"))
    return digest.hexdigest()


def architecture_artifact_sha256(value: Mapping[str, Any]) -> str:
    identity = {
        "source_git_commit": value["source_git_commit"],
        "arm_roster": value["arm_roster"],
        "arms": [
            {
                "arm": arm["arm"],
                "codec": arm["codec"],
                "binary_sha256": arm["binary_sha256"],
                "arm_descriptor_sha256": arm["arm_descriptor_sha256"],
            }
            for arm in value["arms"]
        ],
    }
    digest = hashlib.sha256()
    digest.update(b"wirehair.wh2.architecture-artifacts.v1\0")
    digest.update(canonical_json(identity).encode("utf-8"))
    return digest.hexdigest()


def frozen_arm_artifacts(value: Mapping[str, Any]) -> Mapping[str, Any]:
    return {
        arm["arm"]: {
            "codec": arm["codec"],
            "binary_sha256": arm["binary_sha256"],
            "arm_descriptor_sha256": arm["arm_descriptor_sha256"],
        }
        for arm in value["arms"]
    }


def selected_architecture_sha256(
        arm: str, artifact: Mapping[str, Any]) -> str:
    identity = {
        "arm": arm,
        "codec": artifact["codec"],
        "arm_descriptor_sha256": artifact["arm_descriptor_sha256"],
    }
    digest = hashlib.sha256()
    digest.update(b"wirehair.wh2.selected-architecture.v1\0")
    digest.update(canonical_json(identity).encode("utf-8"))
    return digest.hexdigest()


def contract_sha256(contract: Mapping[str, Any]) -> str:
    """Portable identity over parsed canonical JSON, independent of checkout EOL."""
    return sha256_json(contract)


def serialized_wh2_profile(K: int, block_bytes: int,
                           seed_attempt: int) -> bytes:
    if (type(K) is not int or not 2 <= K <= 64000 or
            type(block_bytes) is not int or block_bytes <= 0 or
            type(seed_attempt) is not int or not 0 <= seed_attempt < 256):
        fail("cannot serialize an invalid WH2 benchmark profile")
    message_bytes = K * block_bytes
    if message_bytes >= 1 << 64 or block_bytes >= 1 << 32:
        fail("WH2 benchmark profile dimensions overflow the public descriptor")
    return b"".join((
        b"WHV2",
        PROFILE_ENCODING_VERSION.to_bytes(2, "little"),
        (32).to_bytes(2, "little"),
        PROFILE_ID.to_bytes(8, "little"),
        message_bytes.to_bytes(8, "little"),
        block_bytes.to_bytes(4, "little"),
        seed_attempt.to_bytes(1, "little"),
        b"\0\0\0",
    ))


def realized_construction_sha256(
        arm_descriptor_sha256: str, K: int, block_bytes: int,
        seed_attempt: int) -> str:
    if (not isinstance(arm_descriptor_sha256, str) or
            SHA256.fullmatch(arm_descriptor_sha256) is None):
        fail("cannot receipt an invalid arm descriptor hash")
    digest = hashlib.sha256()
    digest.update(b"wirehair.wh2.realized-construction.v1\0")
    digest.update(bytes.fromhex(arm_descriptor_sha256))
    digest.update(serialized_wh2_profile(K, block_bytes, seed_attempt))
    return digest.hexdigest()


def _effective_raw_precode_seed(attempt: int) -> str:
    if type(attempt) is not int or not 0 <= attempt < 256:
        fail("raw precode attempt is outside [0,255]")
    seed = (RAW_PRECODE_BASE_SEED +
            attempt * RAW_PRECODE_ATTEMPT_STRIDE) & ((1 << 64) - 1)
    return "0x{:016x}".format(seed)


def _effective_raw_packet_seed(attempt: int) -> str:
    if type(attempt) is not int or not 0 <= attempt < 256:
        fail("raw packet attempt is outside [0,255]")
    seed = (RAW_PACKET_BASE_SEED +
            attempt * RAW_PACKET_ATTEMPT_STRIDE) & ((1 << 32) - 1)
    return "0x{:08x}".format(seed)


def _validate_raw_construction_fields(
        value: Mapping[str, Any], context: str,
        require_dense_anchor_layout: bool = True) -> Mapping[str, Any]:
    fields = RAW_CONSTRUCTION_FIELDS if require_dense_anchor_layout else \
        LEGACY_RAW_CONSTRUCTION_FIELDS
    _exact_keys(value, fields, context)
    if value["construction_seed_basis"] != RAW_CONSTRUCTION_SEED_BASIS:
        fail("{} uses an unknown construction seed basis".format(context))
    if value["seed_schedule_sha256"] != RAW_SEED_SCHEDULE_SHA256:
        fail("{} does not bind the uniform raw seed schedule".format(context))
    precode_attempt = value["precode_attempt"]
    packet_attempt = value["packet_attempt"]
    if (type(precode_attempt) is not int or
            not 0 <= precode_attempt < 256 or
            type(packet_attempt) is not int or
            not 0 <= packet_attempt < 256):
        fail("{} raw attempts must be uint8 integers".format(context))
    precode_seed = value["effective_precode_seed"]
    packet_seed = value["effective_packet_seed"]
    if (not isinstance(precode_seed, str) or
            HEX64.fullmatch(precode_seed) is None or
            precode_seed != _effective_raw_precode_seed(precode_attempt)):
        fail("{} effective precode seed is invalid".format(context))
    if (not isinstance(packet_seed, str) or
            HEX32.fullmatch(packet_seed) is None or
            packet_seed != _effective_raw_packet_seed(packet_attempt)):
        fail("{} effective packet seed is invalid".format(context))
    for field in (
            "staircase", "binary_dense_rows", "gf256_heavy_rows",
            "mix_count"):
        item = value[field]
        if type(item) is not int or not 1 <= item < (1 << 32):
            fail("{} {} must be a positive uint32".format(context, field))
    source_hits = value["source_hits"]
    if type(source_hits) is not int or not 1 <= source_hits <= 8:
        fail("{} source_hits must be an integer in [1,8]".format(context))
    if type(value["dense_identity_corner"]) is not bool:
        fail("{} dense_identity_corner must be boolean".format(context))
    if value["heavy_family"] != "periodic-cauchy":
        fail("{} heavy_family must be periodic-cauchy".format(context))
    if require_dense_anchor_layout and value["dense_anchor_layout"] not in (
            "disabled", "two07"):
        fail("{} dense_anchor_layout must be disabled or two07".format(
            context))
    return value


def _validate_raw_descriptor_structure(
        arm: str, descriptor_sha256: str, K: int,
        value: Mapping[str, Any], context: str) -> None:
    expected = RAW_STRUCTURE_BY_DESCRIPTOR.get(descriptor_sha256)
    if expected is None:
        fail("{} uses an unknown or retired raw descriptor".format(context))
    if arm != expected["arm"]:
        fail("{} raw arm name disagrees with its descriptor".format(context))
    for field in (
            "binary_dense_rows", "gf256_heavy_rows", "heavy_family",
            "mix_count"):
        if value[field] != expected[field]:
            fail("{} {} disagrees with its raw descriptor".format(
                context, field))
    expected_staircase = _raw_staircase_for_K(K)
    if value["staircase"] != expected_staircase:
        fail("{} staircase disagrees with the frozen K geometry".format(
            context))
    if value["source_hits"] != (3 if K >= 10000 else 2):
        fail("{} source_hits disagrees with the frozen K geometry".format(
            context))
    if value["dense_identity_corner"] is not False:
        fail("{} dense identity corner disagrees with the raw descriptor".format(
            context))
    expected_layout = expected["dense_anchor_layout"]
    actual_layout = value.get("dense_anchor_layout", "disabled")
    if actual_layout != expected_layout:
        fail("{} dense anchor layout disagrees with the raw descriptor".format(
            context))


def raw_realized_construction_sha256(
        codec: str, arm: str, arm_descriptor_sha256: str, K: int,
        block_bytes: int, raw_fields: Mapping[str, Any]) -> str:
    if codec != "wirehair2_experiment":
        fail("raw realized construction requires an experimental WH2 codec")
    if (not isinstance(arm, str) or not arm or
            not isinstance(arm_descriptor_sha256, str) or
            SHA256.fullmatch(arm_descriptor_sha256) is None):
        fail("raw realized construction has an invalid arm identity")
    if (type(K) is not int or not 2 <= K <= 64000 or
            type(block_bytes) is not int or block_bytes <= 0 or
            block_bytes >= 1 << 32 or K * block_bytes >= 1 << 64):
        fail("raw realized construction has invalid dimensions")
    fields = _validate_raw_construction_fields(
        raw_fields, "raw realized construction", True)
    _validate_raw_descriptor_structure(
        arm, arm_descriptor_sha256, K, fields,
        "raw realized construction")
    return sha256_json({
        "schema": "wirehair.wh2.raw-realized-construction.v2",
        "codec": codec,
        "arm": arm,
        "arm_descriptor_sha256": arm_descriptor_sha256,
        "K": K,
        "block_bytes": block_bytes,
        **{field: fields[field] for field in RAW_CONSTRUCTION_FIELDS},
    })


def legacy_raw_realized_construction_sha256(
        codec: str, arm: str, arm_descriptor_sha256: str, K: int,
        block_bytes: int, raw_fields: Mapping[str, Any]) -> str:
    """Recompute legacy raw-v2 evidence without permitting Two07 downgrade."""
    if codec != "wirehair2_experiment":
        fail("legacy raw realized construction requires experimental WH2")
    if (not isinstance(arm, str) or not arm.startswith("wirehair2_raw_") or
            not isinstance(arm_descriptor_sha256, str) or
            SHA256.fullmatch(arm_descriptor_sha256) is None):
        fail("legacy raw realized construction has an invalid arm identity")
    structure = RAW_STRUCTURE_BY_DESCRIPTOR.get(arm_descriptor_sha256)
    if structure is None or structure["dense_anchor_layout"] != "disabled":
        fail("legacy raw evidence cannot authenticate a dense-anchor layout")
    if (type(K) is not int or not 2 <= K <= 64000 or
            type(block_bytes) is not int or block_bytes <= 0 or
            block_bytes >= 1 << 32 or K * block_bytes >= 1 << 64):
        fail("legacy raw realized construction has invalid dimensions")
    fields = _validate_raw_construction_fields(
        raw_fields, "legacy raw realized construction", False)
    _validate_raw_descriptor_structure(
        arm, arm_descriptor_sha256, K, fields,
        "legacy raw realized construction")
    return sha256_json({
        "schema": "wirehair.wh2.raw-realized-construction.v1",
        "codec": codec,
        "arm": arm,
        "arm_descriptor_sha256": arm_descriptor_sha256,
        "K": K,
        "block_bytes": block_bytes,
        **{field: fields[field]
           for field in LEGACY_RAW_CONSTRUCTION_FIELDS},
    })


def _load_json_bytes(data: bytes, context: str) -> Any:
    def parse_int(token: str) -> int:
        value = int(token)
        if not -(1 << 63) <= value < (1 << 63):
            fail("{}: JSON integer is outside signed int64".format(context))
        return value

    def parse_float(token: str) -> float:
        value = float(token)
        if not math.isfinite(value):
            fail("{}: non-finite JSON number {}".format(context, token))
        return value

    try:
        text = data.decode("utf-8")
        return json.loads(
            text, object_pairs_hook=_object_no_duplicates,
            parse_int=parse_int, parse_float=parse_float,
            parse_constant=lambda token: fail(
                "{}: nonstandard JSON constant {}".format(context, token)),
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError,
            OverflowError, RecursionError) as exc:
        fail("{}: invalid JSON: {}".format(context, exc))
    return None


def _exact_keys(value: Any, keys: Iterable[str], context: str) -> Mapping[str, Any]:
    expected = set(keys)
    if not isinstance(value, dict) or set(value) != expected:
        fail("{}: expected keys {}, got {}".format(
            context, sorted(expected),
            sorted(value) if isinstance(value, dict) else type(value).__name__,
        ))
    return value


def _integer(value: Any, context: str, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        fail("{}: expected integer >= {}".format(context, minimum))
    return value


def _exact_integer(value: Any, expected: int, context: str) -> int:
    if type(value) is not int or value != expected:
        fail("{}: expected exact integer {}".format(context, expected))
    return value


def _hex_seed(value: Any, context: str) -> str:
    if not isinstance(value, str) or HEX64.fullmatch(value) is None:
        fail("{}: expected canonical 64-bit lowercase hex seed".format(context))
    return value


def _seed_pairs(contract: Mapping[str, Any], mode: str) -> Sequence[Dict[str, Any]]:
    seeds = contract["seeds"]
    if mode == "raw_paired_training":
        base_attempts = seeds["raw_base_seed_attempts"]
        loss_roots = seeds["training_loss_roots"]
    elif mode == "production_training":
        base_attempts = [seeds["production_base_seed_attempt"]] * 3
        loss_roots = seeds["training_loss_roots"]
    elif mode == "production_validation":
        base_attempts = [seeds["production_base_seed_attempt"]] * 3
        loss_roots = seeds["validation_loss_roots"]
    else:
        fail("unknown seed mode: " + str(mode))
    if len(base_attempts) != 3 or len(loss_roots) != 3:
        fail("seed mode does not resolve to exactly three trials")
    return [
        {"base_seed_attempt": base_attempt, "loss_seed": loss_seed}
        for base_attempt, loss_seed in zip(base_attempts, loss_roots)
    ]


def _k_values(contract: Mapping[str, Any], name: str) -> Sequence[int]:
    k_set = contract["k_sets"].get(name)
    if isinstance(k_set, list):
        return k_set
    if isinstance(k_set, dict) and set(k_set) == {"first", "last"}:
        return range(k_set["first"], k_set["last"] + 1)
    fail("unknown or malformed K set: " + name)
    return ()


def _band_for(contract: Mapping[str, Any], K: int) -> str:
    matches = [
        band["name"] for band in contract["k_bands"]
        if band["first"] <= K <= band["last"]
    ]
    if len(matches) != 1:
        fail("K={} is not covered by exactly one band".format(K))
    return matches[0]


def _validate_structure(contract: Any, check_domain_hashes: bool) -> Mapping[str, Any]:
    top = _exact_keys(contract, (
        "schema", "contract_id", "field", "goal", "k_bands", "k_sets",
        "seeds", "strata_sets", "recovery", "timing", "selection",
        "evidence",
    ), "contract")
    if top["schema"] != SCHEMA:
        fail("unexpected contract schema")
    if (top["contract_id"] != "wh2-pure-gf256-1pct-v7" or
            top["field"] != "GF(256)"):
        fail("v4 is pure GF(256) only")

    goal = _exact_keys(top["goal"], (
        "primary_failure_threshold_ppm", "primary_overhead",
        "final_bad_seed_count", "primary_speed_scope",
        "cross_codec_speed_scope",
    ), "goal")
    if (_exact_integer(goal["primary_failure_threshold_ppm"], 10000,
                       "primary failure threshold") != 10000 or
            _exact_integer(goal["primary_overhead"], 0,
                           "primary overhead") != 0 or
            _exact_integer(goal["final_bad_seed_count"], 0,
                           "final bad-seed count") != 0 or
            goal["primary_speed_scope"] != "decoder_solve" or
            goal["cross_codec_speed_scope"] != "receive_to_success"):
        fail("v4 goal constants changed without a schema bump")

    bands = top["k_bands"]
    if not isinstance(bands, list) or len(bands) != 6:
        fail("the contract needs exactly six K bands")
    next_k = 2
    names = set()
    for index, band_value in enumerate(bands):
        band = _exact_keys(band_value, ("name", "first", "last"),
                           "k_bands[{}]".format(index))
        first = _integer(band["first"], "band first", 2)
        last = _integer(band["last"], "band last", first)
        if not isinstance(band["name"], str) or not band["name"] or first != next_k:
            fail("K bands must be named, ordered, and contiguous")
        if band["name"] in names:
            fail("duplicate K band name")
        names.add(band["name"])
        next_k = last + 1
    if next_k != 64001:
        fail("K bands must cover exactly K=2..64000")
    actual_bands = tuple(
        (band["name"], band["first"], band["last"]) for band in bands)
    if actual_bands != EXPECTED_BANDS:
        fail("v4 K bands changed without a schema bump")

    k_sets = _exact_keys(top["k_sets"], ("short", "timing_short", "all"), "k_sets")
    all_set = _exact_keys(k_sets["all"], ("first", "last"), "k_sets.all")
    if (_exact_integer(all_set["first"], 2, "k_sets.all.first") != 2 or
            _exact_integer(all_set["last"], 64000,
                           "k_sets.all.last") != 64000):
        fail("all-K set must be exactly K=2..64000")
    for name in ("short", "timing_short"):
        values = k_sets[name]
        if (not isinstance(values, list) or not values or
                any(type(K) is not int for K in values) or
                values != sorted(set(values)) or values[0] < 2 or values[-1] > 64000):
            fail("{} K set must be a sorted unique in-range list".format(name))
        if {_band_for(top, K) for K in values} != names:
            fail("{} K set must cover all six bands".format(name))
    if tuple(k_sets["short"]) != EXPECTED_SHORT_K or \
            tuple(k_sets["timing_short"]) != EXPECTED_TIMING_SHORT_K:
        fail("v4 bounded K cohorts changed without a schema bump")

    seeds = _exact_keys(top["seeds"], (
        "raw_base_seed_attempts", "production_base_seed_attempt",
        "training_loss_roots", "validation_loss_roots",
    ), "seeds")
    raw_base_attempts = seeds["raw_base_seed_attempts"]
    if (not isinstance(raw_base_attempts, list) or
            len(raw_base_attempts) != 3 or
            any(type(attempt) is not int or not 0 <= attempt < 256
                for attempt in raw_base_attempts) or
            len(set(raw_base_attempts)) != 3):
        fail("raw_base_seed_attempts must contain three unique uint8 values")
    if tuple(raw_base_attempts) != EXPECTED_RAW_BASE_ATTEMPTS:
        fail("v4 raw base seed attempts changed without a schema bump")
    _exact_integer(seeds["production_base_seed_attempt"], 0,
                   "production base seed attempt")
    if seeds["production_base_seed_attempt"] not in raw_base_attempts:
        fail("production base seed attempt must be one raw seed attempt")
    for field in ("training_loss_roots", "validation_loss_roots"):
        roots = seeds[field]
        if not isinstance(roots, list) or len(roots) != 3:
            fail("{} must contain exactly three unique roots".format(field))
        for index, root in enumerate(roots):
            _hex_seed(root, "{}[{}]".format(field, index))
        if len(set(roots)) != 3:
            fail("{} must contain exactly three unique roots".format(field))
    if (tuple(seeds["training_loss_roots"]) != EXPECTED_TRAINING_LOSS_ROOTS or
            tuple(seeds["validation_loss_roots"]) !=
            EXPECTED_VALIDATION_LOSS_ROOTS):
        fail("v4 loss roots changed without a schema bump")
    if not set(seeds["training_loss_roots"]).isdisjoint(
            seeds["validation_loss_roots"]):
        fail("training and validation loss roots must be disjoint")

    strata_sets = _exact_keys(top["strata_sets"], ("development", "hard"), "strata_sets")
    for set_name, expected_count in (("development", 4), ("hard", 3)):
        strata = strata_sets[set_name]
        if not isinstance(strata, list) or len(strata) != expected_count:
            fail("{} needs exactly {} strata".format(set_name, expected_count))
        seen = set()
        for index, stratum_value in enumerate(strata):
            stratum = _exact_keys(stratum_value, ("schedule", "loss_ppm"),
                                  "stratum {}:{}".format(set_name, index))
            key = (stratum["schedule"], _integer(stratum["loss_ppm"], "loss_ppm", 1))
            if (not isinstance(key[0], str) or key[1] >= 1000000 or key in seen):
                fail("loss strata must be unique and have loss in (0,1)")
            seen.add(key)
    actual_strata = {
        set_name: tuple(
            (value["schedule"], value["loss_ppm"])
            for value in strata_sets[set_name])
        for set_name in strata_sets
    }
    if actual_strata != EXPECTED_STRATA_SETS:
        fail("v4 loss/schedule strata changed without a schema bump")

    recovery = _exact_keys(top["recovery"], (
        "overhead_thresholds", "overhead_cap", "raw_construction_attempts",
        "production_seed_fixups_in_raw_phase", "max_construction_attempts",
        "phase_seed_policy", "repair_rule", "attempt_derivation",
        "packet_trace", "cell_key",
        "outcomes", "domains",
    ), "recovery")
    if (recovery["overhead_thresholds"] != [0, 1, 2, 4] or
            any(type(value) is not int
                for value in recovery["overhead_thresholds"]) or
            _exact_integer(recovery["overhead_cap"], 4,
                           "recovery overhead cap") != 4 or
            _exact_integer(recovery["raw_construction_attempts"], 1,
                           "raw construction attempts") != 1 or
            _exact_integer(recovery["production_seed_fixups_in_raw_phase"], 0,
                           "raw production seed fixups") != 0 or
            _exact_integer(recovery["max_construction_attempts"], 256,
                           "maximum construction attempts") != 256):
        fail("v4 raw/overhead constants changed without a schema bump")
    if (recovery["phase_seed_policy"] != EXPECTED_PHASE_SEED_POLICY or
            recovery["repair_rule"] != EXPECTED_REPAIR_RULE or
            recovery["attempt_derivation"] != EXPECTED_ATTEMPT_DERIVATION):
        fail("v4 raw/repaired seed policy changed without a schema bump")
    if recovery["packet_trace"] != EXPECTED_PACKET_TRACE:
        fail("v4 packet-trace algorithm changed without a schema bump")
    expected_key = [
        "phase", "band", "K", "block_bytes", "loss_ppm", "schedule",
        "trial", "base_seed_attempt", "loss_seed", "overhead_cap",
    ]
    if recovery["cell_key"] != expected_key:
        fail("recovery cell key changed without a schema bump")
    outcomes = _exact_keys(recovery["outcomes"], (
        "success", "need_more_at_cap", "construct_failed", "unsupported",
        "fatal_or_internal",
    ), "recovery outcomes")
    if outcomes != EXPECTED_OUTCOMES:
        fail("v4 recovery outcomes changed without a schema bump")

    domains = recovery["domains"]
    expected_phases = {
        "development", "final_raw", "final_repaired", "final_validation",
        "cross_width_validation",
    }
    if not isinstance(domains, dict) or set(domains) != expected_phases:
        fail("recovery phase roster changed without a schema bump")
    for phase in sorted(domains):
        domain = _exact_keys(domains[phase], (
            "k_set", "block_bytes", "seed_mode", "strata_set",
            "expected_cells_per_arm", "domain_sha256",
        ), "recovery domain " + phase)
        expected_k_set, expected_widths, expected_seed_mode, \
            expected_strata_set, expected_count = EXPECTED_RECOVERY_DOMAINS[phase]
        actual_widths = tuple(domain["block_bytes"]) \
            if isinstance(domain["block_bytes"], list) else ()
        if (domain["k_set"] != expected_k_set or
                actual_widths != expected_widths):
            fail("v4 {} K set or widths changed without a schema bump".format(
                phase))
        if (domain["seed_mode"] != expected_seed_mode or
                domain["strata_set"] != expected_strata_set or
                type(domain["expected_cells_per_arm"]) is not int or
                domain["expected_cells_per_arm"] != expected_count):
            fail("v4 {} recovery domain changed without a schema bump".format(
                phase))
        widths = domain["block_bytes"]
        if (not isinstance(widths, list) or not widths or
                any(type(width) is not int or width <= 0 for width in widths) or
                widths != sorted(set(widths))):
            fail("{} block widths must be sorted unique positive integers".format(phase))
        count = (len(_k_values(top, domain["k_set"])) * len(widths) *
                 len(_seed_pairs(top, domain["seed_mode"])) *
                 len(strata_sets[domain["strata_set"]]))
        if domain["expected_cells_per_arm"] != count:
            fail("{} expected cell count is {}, want {}".format(
                phase, domain["expected_cells_per_arm"], count))
        if check_domain_hashes:
            digest = recovery_domain_sha256(top, phase)
            if domain["domain_sha256"] != digest:
                fail("{} domain hash mismatch: want {}".format(phase, digest))

    timing = top["timing"]
    required_timing = {
        "development_wall_time_seconds_per_candidate",
        "loss_ppm", "schedule", "fixed_received_overhead",
        "receive_overhead_cap", "base_cell_key", "cell_key",
        "seed_derivation", "source_derivation", "domains", "confidence",
        "practical_margin_ppm",
        "effective_floor", "aa_repeatability_rule", "order", "cache_state",
        "production_policy", "eligibility", "scope_protocol",
        "required_panels", "scopes", "panel_protocol", "execution_geometry",
        "statistics", "qualification",
        "invalid_timer_policy",
    }
    _exact_keys(timing, required_timing, "timing")
    exact_timing_integer_fields = (
        "development_wall_time_seconds_per_candidate", "loss_ppm",
        "fixed_received_overhead", "receive_overhead_cap",
        "practical_margin_ppm",
    )
    for field in exact_timing_integer_fields:
        _exact_integer(timing[field], EXPECTED_TIMING_POLICY[field],
                       "timing." + field)
    for field in (
            "schedule", "seed_derivation", "source_derivation", "confidence", "effective_floor",
            "aa_repeatability_rule", "order", "cache_state",
            "production_policy", "eligibility", "scope_protocol"):
        if timing[field] != EXPECTED_TIMING_POLICY[field]:
            fail("timing-v6 {} changed without a contract revision".format(
                field))
    if (timing["panel_protocol"] != EXPECTED_PANEL_PROTOCOL or
            timing["execution_geometry"] !=
                EXPECTED_TIMING_EXECUTION_GEOMETRY or
            timing["qualification"] != EXPECTED_TIMING_QUALIFICATION or
            timing["statistics"] != EXPECTED_TIMING_STATISTICS or
            timing["invalid_timer_policy"] != EXPECTED_INVALID_TIMER_POLICY or
            type(timing["panel_protocol"].get("slots_per_panel")) is not int or
            type(timing["panel_protocol"].get(
                "warmups_per_logical_side")) is not int or
            type(timing["panel_protocol"].get(
                "measured_batch_block_target")) is not int or
            timing["statistics"].get("strict_inequalities") is not True or
            any(type(value) is not float for value in
                timing["statistics"].get(
                    "t_critical_by_independent_rounds", {}).values())):
        fail("timing-v6 execution/panel/statistics/qualification protocol changed without "
             "a contract revision")
    if (timing["required_panels"] != [
                "each_arm_AA",
                "candidate_vs_wirehair2_head_decoder_solve",
                "candidate_vs_wirehair1_receive_to_success",
                "candidate_vs_both_controls_encoder"] or
            timing["scopes"] != [
                "decoder_solve", "encoder_init_plus_first_K_symbols",
                "receive_to_success"]):
        fail("timing-v6 constants changed without a contract revision")
    expected_timing_base_key = [
        "phase", "band", "K", "block_bytes", "loss_ppm", "schedule",
        "replicate", "base_seed_attempt", "base_loss_seed",
        "fixed_received_overhead", "receive_overhead_cap",
        "invocations_per_slot", "interleave_policy",
    ]
    expected_timing_key = [
        "phase", "band", "K", "block_bytes", "loss_ppm", "schedule",
        "replicate", "base_seed_attempt", "base_loss_seed",
        "base_cell_sha256", "loss_retry_offset", "loss_seed",
        "fixed_received_overhead", "receive_overhead_cap",
        "invocations_per_slot",
        "interleave_policy",
    ]
    if (timing["base_cell_key"] != expected_timing_base_key or
            timing["cell_key"] != expected_timing_key):
        fail("timing-v6 cell or seed derivation changed without a contract "
             "revision")
    timing_domains = timing["domains"]
    if not isinstance(timing_domains, dict) or set(timing_domains) != {
            "development", "final"}:
        fail("timing domain roster changed without a schema bump")
    expected_timing_domains = {
        "development": (
            "timing_short", [64, 1280], "production_training",
            96, 12, 2304),
        "final": (
            "short", [2, 64, 256, 1280, 4096],
            "production_validation", 96, 12, 14400),
    }
    for phase, expected in expected_timing_domains.items():
        domain = _exact_keys(timing_domains[phase], (
            "k_set", "block_bytes", "seed_mode", "paired_repetitions",
            "independent_rounds", "expected_cells", "base_domain_sha256",
        ), "timing domain " + phase)
        (k_name, expected_widths, seed_mode, repetitions,
         independent_rounds, expected_count) = expected
        if (domain["k_set"] != k_name or domain["block_bytes"] != expected_widths or
                domain["seed_mode"] != seed_mode or
                type(domain["paired_repetitions"]) is not int or
                domain["paired_repetitions"] != repetitions or
                type(domain["independent_rounds"]) is not int or
                domain["independent_rounds"] != independent_rounds or
                type(domain["expected_cells"]) is not int or
                domain["expected_cells"] != expected_count):
            fail("timing-v6 {} domain changed without a contract revision".
                 format(phase))
        widths = domain["block_bytes"]
        if (not isinstance(widths, list) or
                any(type(width) is not int or width <= 0 for width in widths) or
                widths != sorted(set(widths))):
            fail("{} timing widths must be sorted unique positive".format(phase))
        count = len(_k_values(top, k_name)) * len(widths) * repetitions
        if count != expected_count:
            fail("{} timing cell count is inconsistent".format(phase))
        worker_count = timing["execution_geometry"]["timing_worker_count"]
        if (repetitions <= 0 or independent_rounds <= 1 or
                repetitions != independent_rounds * worker_count):
            fail("{} timing repetitions must be exactly eight lanes in each "
                 "independent round".format(phase))
        if check_domain_hashes:
            digest = timing_base_domain_sha256(top, phase)
            if domain["base_domain_sha256"] != digest:
                fail("{} timing base-domain hash mismatch: want {}".format(
                    phase, digest))

    selection = _exact_keys(top["selection"], (
        "controls", "development_architecture_roles",
        "development_recovery_eligibility",
        "development_timing_eligibility",
        "development_timing_proxy_policy", "noninferiority_margin_ppm",
        "architecture_failure_equivalence_ppm", "raw_weak_seed_definition",
        "raw_individual_weak_seeds_are_vetoes",
        "raw_repairs_and_introductions_are_descriptive", "architecture_order",
        "development_resolved_slowdown_rejects",
        "final_decoder_solve_vs_wirehair2",
        "final_receive_to_success_vs_wirehair1", "final_encoder_rule",
        "final_raw_failure_rule", "final_repaired_rule",
        "final_validation_rule",
    ), "selection")
    if (selection != EXPECTED_SELECTION_POLICY or
            type(selection["noninferiority_margin_ppm"]) is not int or
            type(selection["architecture_failure_equivalence_ppm"]) is not int or
            selection["raw_individual_weak_seeds_are_vetoes"] is not False or
            selection["raw_repairs_and_introductions_are_descriptive"] is not True or
            selection["development_resolved_slowdown_rejects"] is not True):
        fail("v4 comparison/weak-seed policy changed without a schema bump")
    evidence = _exact_keys(top["evidence"], (
        "freeze_before_results", "freeze_manifest_schema",
        "trace_manifest_schema", "repair_map_schema",
        "timing_qualification_schema", "timing_qualification_binding",
        "canonical_hashing",
        "common_cell_policy",
        "invalid_domain_mutations", "selectable_excluded_cells",
        "unsupported_policy", "thermal_policy", "final_continuity_policy",
    ), "evidence")
    if (evidence != EXPECTED_EVIDENCE_POLICY or
            type(evidence["selectable_excluded_cells"]) is not int):
        fail("evidence must freeze the roster and permit no excluded cells")
    return top


def load_contract(path: Path = DEFAULT_CONTRACT,
                  check_domain_hashes: bool = True) -> Mapping[str, Any]:
    try:
        data = path.read_bytes()
    except OSError as exc:
        fail("cannot read contract {}: {}".format(path, exc))
    return _validate_structure(_load_json_bytes(data, str(path)), check_domain_hashes)


def _load_canonical_json_file(path: Path, context: str) -> Mapping[str, Any]:
    try:
        data = path.read_bytes()
    except OSError as exc:
        fail("cannot read {} {}: {}".format(context, path, exc))
    value = _load_json_bytes(data, context)
    if not isinstance(value, dict):
        fail("{} must be a JSON object".format(context))
    logical = canonical_json(value).encode("utf-8")
    if data not in (logical + b"\n", logical + b"\r\n"):
        fail("{} must be canonical JSON followed by one line ending".format(
            context))
    return value


def load_freeze_manifest(
        contract: Mapping[str, Any], phase: str, path: Path,
        evidence_kind: str = "recovery",
        timing_qualification: Optional[TimingQualification] = None,
        ) -> Mapping[str, Any]:
    manifest = _load_canonical_json_file(path, "freeze manifest")
    schema = manifest.get("schema")
    manifest_fields = RAW_V3_FREEZE_FIELDS if schema == RAW_FREEZE_SCHEMA \
        else FREEZE_FIELDS
    _exact_keys(manifest, manifest_fields, "freeze manifest")
    if (schema not in (
            FREEZE_SCHEMA, LEGACY_RAW_FREEZE_SCHEMA, RAW_FREEZE_SCHEMA) or
            manifest["contract_sha256"] != contract_sha256(contract) or
            manifest["evidence_kind"] != evidence_kind or
            manifest["phase"] != phase):
        fail("freeze manifest does not bind this contract, kind, and phase")
    domains = contract[evidence_kind]["domains"]
    domain = domains.get(phase)
    if evidence_kind == "timing":
        if timing_qualification is None:
            fail("timing freeze requires a validated qualification map")
        qualification = _require_timing_qualification(
            contract, phase, timing_qualification)
        expected_domain_sha256 = qualification.qualified_domain_sha256
    else:
        qualification = None
        expected_domain_sha256 = domain.get("domain_sha256") \
            if isinstance(domain, Mapping) else None
    if domain is None or manifest["domain_sha256"] != expected_domain_sha256:
        fail("freeze manifest does not bind the frozen domain")
    if (not isinstance(manifest["source_git_commit"], str) or
            GIT_COMMIT.fullmatch(manifest["source_git_commit"]) is None):
        fail("freeze manifest source_git_commit is not a full Git object ID")
    for field in (
            "arm_roster_sha256", "trace_manifest_sha256",
            "repair_training_trace_manifest_sha256"):
        if (not isinstance(manifest[field], str) or
                SHA256.fullmatch(manifest[field]) is None):
            fail("freeze manifest {} is not a SHA-256".format(field))
    roster = manifest["arm_roster"]
    if (not isinstance(roster, list) or
            any(not isinstance(arm, str) or not arm for arm in roster) or
            len(set(roster)) != len(roster) or
            manifest["arm_roster_sha256"] != arm_roster_sha256(roster)):
        fail("freeze manifest has an invalid or unbound arm roster")
    controls = tuple(contract["selection"]["controls"])
    if any(control not in roster for control in controls):
        fail("freeze manifest must contain both declared controls")
    if not any(arm not in controls for arm in roster):
        fail("freeze manifest must contain at least one candidate")
    commands = manifest["commands"]
    if (not isinstance(commands, list) or not commands or
            any(not isinstance(command, list) or not command or
                any(not isinstance(argument, str) or not argument
                    for argument in command)
                for command in commands)):
        fail("freeze manifest commands must be nonempty argv arrays")
    affinity = manifest["cpu_affinity"]
    if (not isinstance(affinity, list) or not affinity or
            any(type(cpu) is not int or cpu < 0 for cpu in affinity) or
            affinity != sorted(set(affinity))):
        fail("freeze manifest CPU affinity must be sorted unique CPU IDs")
    if (not isinstance(manifest["host_identity"], dict) or
            not manifest["host_identity"]):
        fail("freeze manifest host identity must be a nonempty object")
    arm_values = manifest["arms"]
    if not isinstance(arm_values, list) or len(arm_values) != len(roster):
        fail("freeze manifest arm records do not match the roster")
    arms: Dict[str, Mapping[str, Any]] = {}
    raw_phase = evidence_kind == "recovery" and \
        contract["recovery"]["phase_seed_policy"].get(phase) == \
        "one_base_attempt_no_retry"
    if schema == FREEZE_SCHEMA and evidence_kind == "recovery" and \
            phase == "development" and \
            contract.get("schema") == SCHEMA:
        fail("v4 development recovery receipts categorically require a raw freeze schema")
    if schema in (LEGACY_RAW_FREEZE_SCHEMA, RAW_FREEZE_SCHEMA) and \
            (evidence_kind != "recovery" or phase != "development"):
        fail("raw freezes are restricted to development recovery evidence")
    # Development timing uses an unrepaired production-profile attempt, but it
    # is not raw architecture evidence and must never accept the raw v2 schema.
    unrepaired_phase = raw_phase or \
        (evidence_kind == "timing" and phase == "development")
    training_trace_sha = manifest["repair_training_trace_manifest_sha256"]
    if unrepaired_phase and training_trace_sha != "0" * 64:
        fail("unrepaired freezes must use the no-training-trace marker")
    if not unrepaired_phase and training_trace_sha == "0" * 64:
        fail("repaired freezes must bind the training trace manifest")
    if phase == "final_repaired" and \
            training_trace_sha != manifest["trace_manifest_sha256"]:
        fail("the repaired training phase must bind its own frozen traces")
    raw_candidate_count = 0
    if schema == RAW_FREEZE_SCHEMA:
        arm_fields = RAW_V3_FREEZE_ARM_FIELDS
    elif schema == LEGACY_RAW_FREEZE_SCHEMA:
        arm_fields = RAW_FREEZE_ARM_FIELDS
    else:
        arm_fields = FREEZE_ARM_FIELDS
    if schema == RAW_FREEZE_SCHEMA:
        roles = _exact_keys(
            manifest["architecture_roles"], RAW_ARCHITECTURE_ROLE_FIELDS,
            "raw architecture roles")
        if roles != EXPECTED_RAW_ARCHITECTURE_ROLES:
            fail("raw v3 freeze does not bind the exact architecture roles")
        for field in (
                "timing_proxy_witness_sha256", "work_rank_summary_sha256",
                "work_rank_result_stream_sha256", "work_rank_domain_sha256"):
            if (not isinstance(manifest[field], str) or
                    SHA256.fullmatch(manifest[field]) is None or
                    manifest[field] == "0" * 64):
                fail("raw v3 freeze {} is not a nonzero SHA-256".format(
                    field))
    for index, arm_value in enumerate(arm_values):
        arm = _exact_keys(arm_value, arm_fields,
                          "freeze arm {}".format(index))
        name = arm["arm"]
        if name != roster[index] or name in arms:
            fail("freeze arm records must follow the unique frozen roster")
        for field in ("binary_sha256", "arm_descriptor_sha256",
                      "repair_map_sha256"):
            if (not isinstance(arm[field], str) or
                    SHA256.fullmatch(arm[field]) is None):
                fail("freeze arm {} is not a SHA-256".format(field))
        if arm["codec"] not in (
                "wirehair2_certified", "wirehair2_experiment", "wirehair1",
                "routed_composite"):
            fail("freeze arm has an unknown codec kind")
        if (name == "wirehair1") != (arm["codec"] == "wirehair1"):
            fail("codec=wirehair1 is reserved for the Wirehair1 control")
        if (name == "wirehair2_head") != \
                (arm["codec"] == "wirehair2_certified"):
            fail("codec=wirehair2_certified is reserved for the WH2 control")
        if name not in controls and arm["codec"] not in (
                "wirehair2_experiment", "routed_composite"):
            fail("candidate codec must be experimental WH2 or routed composite")
        raw_named = name.startswith("wirehair2_raw_")
        raw_descriptor = arm["arm_descriptor_sha256"] in \
            RAW_ARM_DESCRIPTOR_SHA256S
        if schema == FREEZE_SCHEMA and evidence_kind == "recovery" and \
                phase == "development" and \
                (raw_named or raw_descriptor):
            fail("legacy v1 freezes cannot authenticate raw WH2 arms")
        if schema in (LEGACY_RAW_FREEZE_SCHEMA, RAW_FREEZE_SCHEMA):
            if name not in controls:
                if arm["codec"] != "wirehair2_experiment":
                    fail("raw candidates must be experimental WH2")
                expected_structure = RAW_STRUCTURE_BY_DESCRIPTOR.get(
                    arm["arm_descriptor_sha256"])
                if expected_structure is None:
                    fail("raw candidates must use a current closed descriptor")
                if name != expected_structure["arm"]:
                    fail("raw candidate name disagrees with its descriptor")
                if (schema == LEGACY_RAW_FREEZE_SCHEMA and
                        (not raw_named or
                         expected_structure["dense_anchor_layout"] !=
                            "disabled")):
                    fail("raw v2 cannot authenticate dense-anchor evidence")
                raw_candidate_count += 1
            if arm["codec"] == "routed_composite":
                fail("raw freezes cannot contain routed composites")
        expected_policy = "raw_base" if unrepaired_phase else "repair_map"
        if arm["codec"] == "wirehair1":
            expected_policy = "not_applicable"
        if arm["construction_policy"] != expected_policy:
            fail("freeze arm construction policy disagrees with phase and codec")
        if unrepaired_phase and arm["repair_map_sha256"] != "0" * 64:
            fail("unrepaired freeze arms must use the no-map marker")
        if (not unrepaired_phase and arm["codec"] != "wirehair1" and
                arm["repair_map_sha256"] == "0" * 64):
            fail("repaired Wirehair2/composite arms must freeze a repair map")
        if arm["codec"] == "wirehair1" and \
                arm["repair_map_sha256"] != "0" * 64:
            fail("Wirehair1 must use the no-map marker")
        if schema in (LEGACY_RAW_FREEZE_SCHEMA, RAW_FREEZE_SCHEMA):
            seed_basis = arm["construction_seed_basis"]
            seed_schedule_sha256 = arm["seed_schedule_sha256"]
            if arm["codec"] == "wirehair2_certified":
                expected_basis = PRODUCTION_CONSTRUCTION_SEED_BASIS
                expected_schedule_sha256 = "0" * 64
            elif arm["codec"] == "wirehair1":
                expected_basis = NOT_APPLICABLE_CONSTRUCTION_SEED_BASIS
                expected_schedule_sha256 = "0" * 64
            else:
                expected_basis = RAW_CONSTRUCTION_SEED_BASIS
                expected_schedule_sha256 = RAW_SEED_SCHEDULE_SHA256
            if seed_basis != expected_basis:
                fail("raw v2 freeze arm uses the wrong construction seed basis")
            if (not isinstance(seed_schedule_sha256, str) or
                    SHA256.fullmatch(seed_schedule_sha256) is None or
                    seed_schedule_sha256 != expected_schedule_sha256):
                fail("raw v2 freeze arm uses the wrong seed schedule")
        if schema == RAW_FREEZE_SCHEMA:
            if arm["codec"] == "wirehair1":
                expected_layout = "not-applicable"
            elif arm["codec"] == "wirehair2_certified":
                expected_layout = "disabled"
            else:
                structure = RAW_STRUCTURE_BY_DESCRIPTOR.get(
                    arm["arm_descriptor_sha256"])
                expected_layout = None if structure is None else \
                    structure["dense_anchor_layout"]
            if arm["dense_anchor_layout"] != expected_layout:
                fail("raw v3 freeze arm uses the wrong dense-anchor layout")
        arms[name] = arm
    if schema in (LEGACY_RAW_FREEZE_SCHEMA, RAW_FREEZE_SCHEMA) and \
            raw_candidate_count == 0:
        fail("raw freeze must contain at least one raw candidate")
    if schema == RAW_FREEZE_SCHEMA:
        expected_roster = roles["descriptive_controls"] + [
            roles["recovery_reference"]] + roles["recovery_candidates"]
        if roster != expected_roster or roles["timing_proxy"] not in \
                roles["descriptive_controls"]:
            fail("raw v3 freeze roster disagrees with its architecture roles")
    if qualification is not None:
        if qualification.source_git_commit != manifest["source_git_commit"]:
            fail("timing qualification substitutes the frozen source commit")
        for control_tuple in qualification.controls:
            control = dict(zip(
                TIMING_QUALIFICATION_CONTROL_ORDER, control_tuple))
            frozen = arms.get(control["arm"])
            if frozen is None or any(
                    frozen[field] != control[field]
                    for field in (
                        "binary_sha256", "arm_descriptor_sha256",
                        "construction_policy", "repair_map_sha256")):
                fail("timing qualification substitutes a frozen control")
    result = dict(manifest)
    result["arms_by_name"] = arms
    return result


def validate_selection_receipt(
        contract: Mapping[str, Any], value: Mapping[str, Any]) -> Mapping[str, Any]:
    _exact_keys(value, SELECTION_FIELDS, "architecture selection receipt")
    unsigned = {key: item for key, item in value.items()
                if key != "selection_sha256"}
    for field in (
            "selection_sha256", "recovery_freeze_manifest_sha256",
            "timing_freeze_manifest_sha256", "architecture_artifact_sha256",
            "timing_domain_sha256", "timing_qualification_map_sha256",
            "selected_arm_descriptor_sha256",
            "selected_architecture_sha256"):
        if (not isinstance(value[field], str) or
                SHA256.fullmatch(value[field]) is None):
            fail("architecture selection {} is not a SHA-256".format(field))
    if (value["schema"] != SCHEMA + ".architecture-selection.v1" or
            value["contract_sha256"] != contract_sha256(contract) or
            value["recovery_domain_sha256"] != contract["recovery"][
                "domains"]["development"]["domain_sha256"] or
            value["timing_base_domain_sha256"] != contract["timing"][
                "domains"]["development"]["base_domain_sha256"] or
            value["selection_sha256"] != sha256_json(unsigned)):
        fail("architecture selection receipt identity is invalid")
    _exact_integer(
        value["recovery_cells_per_arm"], contract["recovery"]["domains"]
        ["development"]["expected_cells_per_arm"],
        "architecture selection recovery cells")
    candidate_roster = value["candidate_roster"]
    controls = set(contract["selection"]["controls"])
    if (not isinstance(candidate_roster, list) or not candidate_roster or
            any(not isinstance(arm, str) or not arm
                for arm in candidate_roster) or
            candidate_roster != sorted(set(candidate_roster)) or
            set(candidate_roster) & controls):
        fail("architecture selection candidate roster is malformed")
    protocol = contract["timing"]["panel_protocol"]
    panel_count = len(protocol["control_aa"]) + len(candidate_roster) * (
        len(protocol["candidate_aa_scopes"]) + len(protocol["candidate_ab"]))
    expected_timing_rows = contract["timing"]["domains"]["development"][
        "expected_cells"] * panel_count
    _exact_integer(value["timing_rows"], expected_timing_rows,
                   "architecture selection timing rows")
    eligible = value["eligible_candidates"]
    eligible_failures = value["eligible_overhead0_failures"]
    equivalent = value["recovery_equivalent_candidates"]
    if (not isinstance(eligible, list) or
            any(not isinstance(arm, str) or not arm for arm in eligible) or
            eligible != sorted(set(eligible)) or
            not isinstance(equivalent, list) or
            any(not isinstance(arm, str) or not arm for arm in equivalent) or
            equivalent != sorted(set(equivalent)) or
            not isinstance(eligible_failures, dict) or
            set(eligible_failures) != set(eligible) or
            any(type(count) is not int for count in eligible_failures.values()) or
            not set(eligible).issubset(candidate_roster) or
            not set(equivalent).issubset(eligible)):
        fail("architecture selection candidate sets are malformed")
    selected = value["selected_arm"]
    selected_codec = value["selected_codec"]
    minimum = value["minimum_overhead0_failures"]
    allowance = value["recovery_equivalence_allowance"]
    ranking = value["ranking"]
    cells = contract["recovery"]["domains"]["development"][
        "expected_cells_per_arm"]
    margin_ppm = contract["selection"][
        "architecture_failure_equivalence_ppm"]
    expected_allowance = max(
        1, (margin_ppm * cells + 1000000 - 1) // 1000000)
    if (not isinstance(selected, str) or not selected or
            selected not in equivalent or type(minimum) is not int or
            not 0 <= minimum <= cells or type(allowance) is not int or
            allowance != expected_allowance or
            selected_codec not in ("wirehair2_experiment", "routed_composite") or
            not isinstance(ranking, list) or not ranking):
        fail("architecture selection has no promotable winner")
    if (any(not 0 <= count <= cells for count in eligible_failures.values()) or
            minimum != min(eligible_failures.values()) or
            set(equivalent) != {
                arm for arm, count in eligible_failures.items()
                if count <= minimum + allowance}):
        fail("architecture recovery-equivalent set is invalid")
    selected_artifact = {
        "codec": selected_codec,
        "arm_descriptor_sha256": value["selected_arm_descriptor_sha256"],
    }
    if value["selected_architecture_sha256"] != \
            selected_architecture_sha256(selected, selected_artifact):
        fail("selected architecture identity is invalid")
    ranking_values = []
    for index, item_value in enumerate(ranking):
        item = _exact_keys(item_value, (
            "arm", "decoder_solve_mean_log_ratio", "failures_overhead0",
            "failures_overhead1", "failures_overhead2", "failures_overhead4",
        ), "architecture ranking item {}".format(index))
        mean = item["decoder_solve_mean_log_ratio"]
        counts = tuple(item[field] for field in (
            "failures_overhead0", "failures_overhead1",
            "failures_overhead2", "failures_overhead4"))
        if (not isinstance(item["arm"], str) or type(mean) is not float or
                not math.isfinite(mean) or
                any(type(count) is not int or not 0 <= count <= cells
                    for count in counts) or
                any(left < right for left, right in zip(counts, counts[1:]))):
            fail("architecture ranking item is malformed")
        if eligible_failures.get(item["arm"]) != counts[0]:
            fail("architecture ranking overhead-zero count is inconsistent")
        ranking_values.append((
            mean, counts[1], counts[2], counts[3], item["arm"]))
    ranked_arms = [item[-1] for item in ranking_values]
    if (ranking_values != sorted(ranking_values) or
            len(set(ranked_arms)) != len(ranked_arms) or
            set(ranked_arms) != set(equivalent) or
            ranking_values[0][-1] != selected):
        fail("architecture ranking disagrees with the selected winner")
    return value


def validate_final_freeze_continuity(
        contract: Mapping[str, Any],
        freeze_paths: Mapping[Tuple[str, str], Path],
        selection_receipt: Mapping[str, Any],
        timing_qualification: TimingQualification,
        timing_trace_manifest_path: Path,
        ) -> Mapping[str, Any]:
    selection = validate_selection_receipt(contract, selection_receipt)
    if set(freeze_paths) != FINAL_FREEZE_KEYS:
        fail("final continuity requires exactly the five frozen final phases")
    freezes = {
        key: load_freeze_manifest(
            contract, key[1], path, key[0],
            timing_qualification if key == ("timing", "final") else None)
        for key, path in freeze_paths.items()
    }
    final_timing_freeze = freezes[("timing", "final")]
    load_timing_trace_manifest(
        contract, "final", timing_trace_manifest_path,
        final_timing_freeze["trace_manifest_sha256"],
        timing_qualification)
    anchor = freezes[("recovery", "final_raw")]
    controls = tuple(contract["selection"]["controls"])
    candidates = tuple(
        arm for arm in anchor["arm_roster"] if arm not in controls)
    if len(candidates) != 1 or len(anchor["arm_roster"]) != len(controls) + 1:
        fail("final phases must contain both controls and one selected candidate")
    if candidates[0] != selection["selected_arm"]:
        fail("final phases substitute the selected development architecture")
    final_selected_artifact = frozen_arm_artifacts(anchor)[candidates[0]]
    if selected_architecture_sha256(
            candidates[0], final_selected_artifact) != \
            selection["selected_architecture_sha256"]:
        fail("final phases substitute the selected architecture descriptor")
    artifact_sha256 = architecture_artifact_sha256(anchor)
    host_identity = canonical_json(anchor["host_identity"])
    for key, freeze in freezes.items():
        if architecture_artifact_sha256(freeze) != artifact_sha256:
            fail("{}:{} substitutes final architecture artifacts".format(*key))
        if canonical_json(freeze["host_identity"]) != host_identity:
            fail("{}:{} substitutes the final benchmark host".format(*key))

    repaired = freezes[("recovery", "final_repaired")]
    training_trace_sha256 = repaired["trace_manifest_sha256"]
    production_keys = FINAL_FREEZE_KEYS - {("recovery", "final_raw")}
    production_maps = {
        arm: repaired["arms_by_name"][arm]["repair_map_sha256"]
        for arm in repaired["arm_roster"]
    }
    for key in production_keys:
        freeze = freezes[key]
        if freeze["repair_training_trace_manifest_sha256"] != \
                training_trace_sha256:
            fail("{}:{} substitutes the repair-training trace".format(*key))
        actual_maps = {
            arm: freeze["arms_by_name"][arm]["repair_map_sha256"]
            for arm in freeze["arm_roster"]
        }
        if actual_maps != production_maps:
            fail("{}:{} substitutes production repair maps".format(*key))

    identity = {
        "contract_sha256": contract_sha256(contract),
        "architecture_artifact_sha256": artifact_sha256,
        "host_identity": anchor["host_identity"],
        "repair_training_trace_manifest_sha256": training_trace_sha256,
        "production_repair_maps": production_maps,
        "selection_sha256": selection["selection_sha256"],
    }
    return {
        "schema": SCHEMA + ".final-continuity-summary.v1",
        "selected_candidate": candidates[0],
        "selection_sha256": selection["selection_sha256"],
        "promotion_identity_sha256": sha256_json(identity),
        "freeze_manifest_sha256": {
            "{}:{}".format(*key): freeze_manifest_sha256(freezes[key])
            for key in sorted(FINAL_FREEZE_KEYS)
        },
    }


def load_repair_map(
        contract: Mapping[str, Any], path: Path,
        freeze_arm: Mapping[str, Any], source_git_commit: str,
        training_trace_manifest_sha256: str) -> Mapping[int, int]:
    value = _load_canonical_json_file(path, "repair map")
    _exact_keys(value, REPAIR_MAP_FIELDS, "repair map")
    if (value["schema"] != REPAIR_MAP_SCHEMA or
            value["contract_sha256"] != contract_sha256(contract) or
            value["training_domain_sha256"] !=
            contract["recovery"]["domains"]["final_repaired"]["domain_sha256"] or
            value["arm"] != freeze_arm["arm"] or
            value["source_git_commit"] != source_git_commit or
            value["binary_sha256"] != freeze_arm["binary_sha256"] or
            value["arm_descriptor_sha256"] !=
            freeze_arm["arm_descriptor_sha256"] or
            value["training_trace_manifest_sha256"] !=
            training_trace_manifest_sha256 or
            value["entry_kind"] != "retry_offset_indexed_by_K_minus_2" or
            value["attempt_derivation"] != EXPECTED_ATTEMPT_DERIVATION or
            value["repair_rule"] != EXPECTED_REPAIR_RULE):
        fail("repair map does not bind the frozen training artifact")
    _exact_integer(
        value["production_base_seed_attempt"],
        contract["seeds"]["production_base_seed_attempt"],
        "repair map production base seed attempt")
    if repair_map_sha256(value) != freeze_arm["repair_map_sha256"]:
        fail("repair map SHA-256 differs from the pre-result freeze")
    retry_offsets = value["retry_offsets"]
    K_values = _k_values(contract, "all")
    if (not isinstance(retry_offsets, list) or
            len(retry_offsets) != len(K_values)):
        fail("repair map must contain one retry offset for every K")
    attempts: Dict[int, int] = {}
    for expected_K, retry_offset in zip(K_values, retry_offsets):
        retry_offset = _integer(
            retry_offset, "repair map retry offset for K={}".format(expected_K))
        if retry_offset >= contract["recovery"]["max_construction_attempts"]:
            fail("repair map retry offset exceeds uint8 range")
        attempts[expected_K] = (
            value["production_base_seed_attempt"] + retry_offset) & 0xff
    return attempts


def generic_realized_construction_sha256(
        codec: str, arm_descriptor_sha256: str, K: int, block_bytes: int,
        construction_attempt: int) -> str:
    if codec == "wirehair2_certified":
        return realized_construction_sha256(
            arm_descriptor_sha256, K, block_bytes, construction_attempt)
    return sha256_json({
        "schema": "wirehair.wh2.realized-construction.v1",
        "codec": codec,
        "arm_descriptor_sha256": arm_descriptor_sha256,
        "K": K,
        "block_bytes": block_bytes,
        "construction_attempt": construction_attempt,
    })


def expected_construction_attempt(
        frozen_arm: Mapping[str, Any], K: int, base_seed_attempt: int,
        repair_attempts: Mapping[str, Mapping[int, int]]) -> int:
    arm = frozen_arm["arm"]
    if arm in repair_attempts:
        return repair_attempts[arm][K]
    if frozen_arm["construction_policy"] == "not_applicable":
        return 0
    if frozen_arm["construction_policy"] == "raw_base":
        return base_seed_attempt
    fail("mapped construction policy has no authenticated repair map")
    return 0


def iter_recovery_cells(contract: Mapping[str, Any], phase: str) -> Iterator[Dict[str, Any]]:
    domains = contract["recovery"]["domains"]
    if phase not in domains:
        fail("unknown recovery phase: " + phase)
    domain = domains[phase]
    K_values = _k_values(contract, domain["k_set"])
    for block_bytes in domain["block_bytes"]:
        for trial, seed_pair in enumerate(
                _seed_pairs(contract, domain["seed_mode"])):
            for stratum in contract["strata_sets"][domain["strata_set"]]:
                for K in K_values:
                    yield {
                        "phase": phase,
                        "band": _band_for(contract, K),
                        "K": K,
                        "block_bytes": block_bytes,
                        "loss_ppm": stratum["loss_ppm"],
                        "schedule": stratum["schedule"],
                        "trial": trial,
                        "base_seed_attempt": seed_pair["base_seed_attempt"],
                        "loss_seed": seed_pair["loss_seed"],
                        "overhead_cap": contract["recovery"]["overhead_cap"],
                    }


def recovery_domain_sha256(contract: Mapping[str, Any], phase: str) -> str:
    digest = hashlib.sha256()
    for cell in iter_recovery_cells(contract, phase):
        digest.update(canonical_json(cell).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _splitmix64(value: int) -> int:
    mask = (1 << 64) - 1
    value = (value + 0x9e3779b97f4a7c15) & mask
    value = ((value ^ (value >> 30)) * 0xbf58476d1ce4e5b9) & mask
    value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & mask
    return (value ^ (value >> 31)) & mask


def timing_invocations_per_slot(
        contract: Mapping[str, Any], K: int) -> int:
    """Return frozen n, split across each pair of opposite-order slots."""
    if type(K) is not int or K <= 0:
        fail("timing batch K must be a positive integer")
    timing = contract.get("timing")
    protocol = timing.get("panel_protocol") \
        if isinstance(timing, Mapping) else None
    target = protocol.get("measured_batch_block_target") \
        if isinstance(protocol, Mapping) else None
    formula = protocol.get("invocations_per_slot") \
        if isinstance(protocol, Mapping) else None
    if (type(target) is not int or
            target != MEASURED_BATCH_BLOCK_TARGET or
            formula != INVOCATIONS_PER_SLOT_FORMULA):
        fail("timing measured-batch protocol is invalid")
    return max(2, (target + K - 1) // K)


def timing_invocations_for_elapsed_slot(
        contract: Mapping[str, Any], K: int, slot: int) -> int:
    """Return the exact fresh-call count for elapsed slot 0 through 7."""
    if type(slot) is not int or not 0 <= slot < 8:
        fail("timing elapsed slot must be an integer in [0,7]")
    invocations = timing_invocations_per_slot(contract, K)
    return (invocations + 1) // 2 if slot < 4 else invocations // 2


def timing_cohort_count(
        contract: Mapping[str, Any], phase: str,
        frozen_arms: Sequence[str]) -> int:
    """Derive the exact phase cohort count for one frozen arm roster."""
    timing = contract.get("timing") if isinstance(contract, Mapping) else None
    geometry = timing.get("execution_geometry") \
        if isinstance(timing, Mapping) else None
    domains = timing.get("domains") if isinstance(timing, Mapping) else None
    if not isinstance(phase, str):
        fail("timing phase must be a string")
    domain = domains.get(phase) if isinstance(domains, Mapping) else None
    if geometry != EXPECTED_TIMING_EXECUTION_GEOMETRY or \
            not isinstance(domain, Mapping):
        fail("timing execution geometry or phase is invalid")
    if (not isinstance(frozen_arms, Sequence) or
            isinstance(frozen_arms, (str, bytes)) or not frozen_arms or
            any(not isinstance(arm, str) or not arm for arm in frozen_arms) or
            len(frozen_arms) != len(set(frozen_arms))):
        fail("timing frozen arm roster must be nonempty unique strings")
    controls = set(contract["selection"]["controls"])
    if not controls.issubset(set(frozen_arms)):
        fail("timing frozen arm roster omits a required control")
    repetitions = domain.get("paired_repetitions")
    independent_rounds = domain.get("independent_rounds")
    expected_cells = domain.get("expected_cells")
    if (type(repetitions) is not int or repetitions <= 0 or
            type(independent_rounds) is not int or independent_rounds <= 1 or
            repetitions != independent_rounds * geometry["jobs_per_wave"] or
            type(expected_cells) is not int or expected_cells <= 0 or
            expected_cells % repetitions != 0):
        fail("timing phase cardinality cannot form exact cohorts")
    panels = timing_panels(contract, frozen_arms)
    if not panels:
        fail("timing frozen arm roster produces no panels")
    return expected_cells // repetitions * len(panels)


def timing_worker_slot(
        contract: Mapping[str, Any], phase: str,
        frozen_arms: Sequence[str], cohort_index: int,
        replicate: int) -> int:
    """Return the frozen sorted-worker slot for one cohort replicate."""
    timing = contract.get("timing") if isinstance(contract, Mapping) else None
    geometry = timing.get("execution_geometry") \
        if isinstance(timing, Mapping) else None
    domains = timing.get("domains") if isinstance(timing, Mapping) else None
    domain = domains.get(phase) if isinstance(domains, Mapping) else None
    cohort_count = timing_cohort_count(contract, phase, frozen_arms)
    if type(cohort_index) is not int or not 0 <= cohort_index < cohort_count:
        fail("timing cohort index is outside the frozen phase/arm roster")
    repetitions = domain.get("paired_repetitions") \
        if isinstance(domain, Mapping) else None
    if (type(replicate) is not int or type(repetitions) is not int or
            not 0 <= replicate < repetitions):
        fail("timing replicate is outside the frozen phase")
    workers = geometry["timing_worker_count"]
    jobs_per_wave = geometry["jobs_per_wave"]
    independent_round = replicate // jobs_per_wave
    lane = replicate % jobs_per_wave
    return (lane + cohort_index + independent_round) % workers


def iter_timing_base_cells(
        contract: Mapping[str, Any], phase: str) -> Iterator[Dict[str, Any]]:
    """Yield the arm- and outcome-free timing domain before loss retry."""
    domains = contract["timing"]["domains"]
    if phase not in domains:
        fail("unknown timing phase: " + phase)
    domain = domains[phase]
    seed_pairs = _seed_pairs(contract, domain["seed_mode"])
    mask = (1 << 64) - 1
    lanes = contract["timing"]["execution_geometry"]["jobs_per_wave"]
    rounds = domain["independent_rounds"]
    if domain["paired_repetitions"] != lanes * rounds:
        fail("timing domain does not contain eight lanes per round")
    for independent_round in range(rounds):
        for block_bytes in domain["block_bytes"]:
            for K in _k_values(contract, domain["k_set"]):
                for lane in range(lanes):
                    replicate = independent_round * lanes + lane
                    pair = seed_pairs[replicate % len(seed_pairs)]
                    salt = (replicate * 0x9e3779b97f4a7c15) & mask
                    base_loss_seed = _splitmix64(
                        int(pair["loss_seed"], 16) ^ salt)
                    yield {
                        "phase": phase,
                        "band": _band_for(contract, K),
                        "K": K,
                        "block_bytes": block_bytes,
                        "loss_ppm": contract["timing"]["loss_ppm"],
                        "schedule": contract["timing"]["schedule"],
                        "replicate": replicate,
                        "base_seed_attempt": pair["base_seed_attempt"],
                        "base_loss_seed":
                            "0x{:016x}".format(base_loss_seed),
                        "fixed_received_overhead":
                            contract["timing"]["fixed_received_overhead"],
                        "receive_overhead_cap":
                            contract["timing"]["receive_overhead_cap"],
                        "invocations_per_slot":
                            timing_invocations_per_slot(contract, K),
                        "interleave_policy": TIMING_INTERLEAVE_POLICY,
                    }


def timing_base_domain_sha256(
        contract: Mapping[str, Any], phase: str) -> str:
    digest = hashlib.sha256()
    for cell in iter_timing_base_cells(contract, phase):
        digest.update(canonical_json(cell).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _qualified_timing_loss_seed(base_loss_seed: str, retry_offset: int) -> str:
    base = int(_hex_seed(base_loss_seed, "timing base loss seed"), 16)
    if (type(retry_offset) is not int or
            not 0 <= retry_offset < TIMING_LOSS_RETRY_ATTEMPTS):
        fail("timing loss retry offset must be a uint8 integer")
    realized = (base + retry_offset * TIMING_LOSS_RETRY_STRIDE) & \
        ((1 << 64) - 1)
    return "0x{:016x}".format(realized)


def _iter_qualified_timing_cells(
        contract: Mapping[str, Any], phase: str,
        retry_offsets: Sequence[int]) -> Iterator[Dict[str, Any]]:
    expected_cells = contract["timing"]["domains"][phase]["expected_cells"]
    if (not isinstance(retry_offsets, Sequence) or
            isinstance(retry_offsets, (str, bytes)) or
            len(retry_offsets) != expected_cells):
        fail("timing qualification must contain one offset per base cell")
    for ordinal, base_cell in enumerate(iter_timing_base_cells(
            contract, phase)):
        retry_offset = retry_offsets[ordinal]
        loss_seed = _qualified_timing_loss_seed(
            base_cell["base_loss_seed"], retry_offset)
        cell = dict(base_cell)
        cell["base_cell_sha256"] = sha256_json(base_cell)
        cell["loss_retry_offset"] = retry_offset
        cell["loss_seed"] = loss_seed
        yield cell


def _timing_domain_sha256_from_offsets(
        contract: Mapping[str, Any], phase: str,
        retry_offsets: Sequence[int]) -> str:
    digest = hashlib.sha256()
    for cell in _iter_qualified_timing_cells(
            contract, phase, retry_offsets):
        digest.update(canonical_json(cell).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _require_timing_qualification(
        contract: Mapping[str, Any], phase: str,
        qualification: TimingQualification) -> TimingQualification:
    if not _timing_qualification_is_registered(qualification):
        fail("timing qualification object did not come from the strict loader "
             "or was mutated after loading")
    domain = contract["timing"]["domains"].get(phase)
    if (not isinstance(domain, Mapping) or
            qualification.contract_sha256 != contract_sha256(contract) or
            qualification.phase != phase or
            qualification.base_domain_sha256 !=
                domain["base_domain_sha256"] or
            len(qualification.retry_offsets) != domain["expected_cells"] or
            len(qualification.selected_trace_sha256s) !=
                domain["expected_cells"] or
            not isinstance(qualification.map_sha256, str) or
            SHA256.fullmatch(qualification.map_sha256) is None or
            not isinstance(qualification.qualified_domain_sha256, str) or
            SHA256.fullmatch(qualification.qualified_domain_sha256) is None or
            not isinstance(qualification.qualification_audit_sha256, str) or
            SHA256.fullmatch(
                qualification.qualification_audit_sha256) is None or
            not isinstance(qualification.selected_trace_roster_sha256, str) or
            SHA256.fullmatch(
                qualification.selected_trace_roster_sha256) is None or
            not isinstance(qualification.source_git_commit, str) or
            GIT_COMMIT.fullmatch(qualification.source_git_commit) is None):
        fail("timing qualification does not bind this phase and base domain")
    if any(not isinstance(trace_sha256, str) or
           SHA256.fullmatch(trace_sha256) is None
           for trace_sha256 in qualification.selected_trace_sha256s):
        fail("timing qualification has a malformed selected trace identity")
    if timing_selected_trace_roster_sha256(
            qualification.selected_trace_sha256s) != \
            qualification.selected_trace_roster_sha256:
        fail("timing qualification object has a forged selected trace roster")
    if (not isinstance(qualification.controls, tuple) or
            any(not isinstance(control, tuple) or
                len(control) != len(TIMING_QUALIFICATION_CONTROL_ORDER) or
                any(not isinstance(item, str) for item in control)
                for control in qualification.controls)):
        fail("timing qualification has a mutable or malformed control roster")
    control_records = [
        dict(zip(TIMING_QUALIFICATION_CONTROL_ORDER, control))
        for control in qualification.controls
    ]
    _validate_timing_qualification_controls(
        contract, phase, control_records)
    if any(type(offset) is not int or
           not 0 <= offset < TIMING_LOSS_RETRY_ATTEMPTS
           for offset in qualification.retry_offsets):
        fail("timing qualification contains a non-uint8 retry offset")
    qualified_digest = _timing_domain_sha256_from_offsets(
        contract, phase, qualification.retry_offsets)
    if qualified_digest != qualification.qualified_domain_sha256:
        fail("timing qualification object has a forged qualified domain")
    canonical_map = {
        "schema": TIMING_QUALIFICATION_MAP_SCHEMA,
        "contract_sha256": qualification.contract_sha256,
        "phase": qualification.phase,
        "source_git_commit": qualification.source_git_commit,
        "base_domain_sha256": qualification.base_domain_sha256,
        "qualified_domain_sha256": qualification.qualified_domain_sha256,
        "entry_kind": TIMING_QUALIFICATION_ENTRY_KIND,
        "controls": control_records,
        "qualification_audit_sha256":
            qualification.qualification_audit_sha256,
        "selected_trace_roster_sha256":
            qualification.selected_trace_roster_sha256,
        "retry_offsets": list(qualification.retry_offsets),
    }
    if timing_qualification_map_sha256(canonical_map) != \
            qualification.map_sha256:
        fail("timing qualification object has a forged map SHA-256")
    return qualification


def iter_timing_cells(
        contract: Mapping[str, Any], phase: str,
        qualification: TimingQualification) -> Iterator[Dict[str, Any]]:
    """Yield the arm-free qualified cells consumed identically by every arm."""
    validated = _require_timing_qualification(contract, phase, qualification)
    yield from _iter_qualified_timing_cells(
        contract, phase, validated.retry_offsets)


def timing_domain_sha256(
        contract: Mapping[str, Any], phase: str,
        qualification: TimingQualification) -> str:
    validated = _require_timing_qualification(contract, phase, qualification)
    return _timing_domain_sha256_from_offsets(
        contract, phase, validated.retry_offsets)


def timing_panels(
        contract: Mapping[str, Any], arms: Sequence[str],
        ) -> Sequence[Dict[str, str]]:
    controls = tuple(contract["selection"]["controls"])
    candidates = tuple(arm for arm in arms if arm not in controls)
    panels: List[Dict[str, str]] = []
    protocol = contract["timing"]["panel_protocol"]
    for value in protocol["control_aa"]:
        panels.append({
            "panel_kind": "AA", "scope": value["scope"],
            "left_arm": value["arm"], "right_arm": value["arm"],
        })
    for candidate in candidates:
        for scope in protocol["candidate_aa_scopes"]:
            panels.append({
                "panel_kind": "AA", "scope": scope,
                "left_arm": candidate, "right_arm": candidate,
            })
        for value in protocol["candidate_ab"]:
            panels.append({
                "panel_kind": "AB", "scope": value["scope"],
                "left_arm": candidate, "right_arm": value["control"],
            })
    if len({canonical_json(panel) for panel in panels}) != len(panels):
        fail("internal duplicate timing panel")
    return panels


def timing_order(panel: Mapping[str, str], replicate: int) -> str:
    phase_bit = bytes.fromhex(sha256_json(panel))[-1] & 1
    return "ABBA" if (replicate & 1) == phase_bit else "BAAB"


def _timing_four_slot_log_contrast(
        logs: Sequence[float], start: int, order: str) -> float:
    """Return one left/right-oriented four-slot paired log contrast."""
    if len(logs) != 8 or start not in (0, 4):
        fail("timing log contrast requires one block of an eight-slot job")
    if order == "ABBA":
        return ((logs[start] - logs[start + 1]) +
                (logs[start + 3] - logs[start + 2])) / 2.0
    if order == "BAAB":
        return ((logs[start + 1] - logs[start]) +
                (logs[start + 2] - logs[start + 3])) / 2.0
    fail("timing log contrast has an unknown order")
    return 0.0


def _timing_cell_indexes(
        contract: Mapping[str, Any], phase: str,
        qualification: TimingQualification,
        ) -> Mapping[str, Tuple[int, Mapping[str, Any]]]:
    result: Dict[str, Tuple[int, Mapping[str, Any]]] = {}
    for ordinal, cell in enumerate(iter_timing_cells(
            contract, phase, qualification)):
        key = canonical_json(cell)
        if key in result:
            fail("internal duplicate timing cell")
        result[key] = (ordinal, cell)
    return result


def _timing_cell_ordinal(
        contract: Mapping[str, Any], phase: str, row: Mapping[str, Any],
        indexes: Mapping[str, Tuple[int, Mapping[str, Any]]],
        ) -> Tuple[int, Mapping[str, Any]]:
    integer_fields = (
        "K", "block_bytes", "loss_ppm", "replicate", "base_seed_attempt",
        "loss_retry_offset", "fixed_received_overhead", "receive_overhead_cap",
        "invocations_per_slot",
    )
    string_fields = (
        "phase", "band", "schedule", "base_loss_seed", "base_cell_sha256",
        "loss_seed", "interleave_policy",
    )
    if (any(type(row[field]) is not int for field in integer_fields) or
            any(not isinstance(row[field], str) for field in string_fields)):
        fail("timing cell key uses a noncanonical scalar type")
    key_value = {
        key: row[key] for key in contract["timing"]["cell_key"]
    }
    try:
        return indexes[canonical_json(key_value)]
    except KeyError:
        fail("timing row is outside the frozen {} domain".format(phase))
    return 0, {}


def _domain_indexes(contract: Mapping[str, Any], phase: str) -> Tuple[
        Mapping[int, int], Mapping[Tuple[str, str], int],
        Mapping[Tuple[str, int], int], Mapping[int, int]]:
    domain = contract["recovery"]["domains"][phase]
    return (
        {K: index for index, K in enumerate(_k_values(contract, domain["k_set"]))},
        {(pair["base_seed_attempt"], pair["loss_seed"]): index
         for index, pair in enumerate(_seed_pairs(contract, domain["seed_mode"]))},
        {(value["schedule"], value["loss_ppm"]): index
         for index, value in enumerate(contract["strata_sets"][domain["strata_set"]])},
        {width: index for index, width in enumerate(domain["block_bytes"])},
    )


def _cell_ordinal(
        contract: Mapping[str, Any], phase: str, row: Mapping[str, Any],
        indexes: Tuple[Mapping[int, int], Mapping[Tuple[str, str], int],
                       Mapping[Tuple[str, int], int], Mapping[int, int]],
) -> Tuple[int, Dict[str, Any]]:
    domain = contract["recovery"]["domains"][phase]
    k_index, seed_index, stratum_index, width_index = indexes
    integer_fields = (
        "K", "block_bytes", "loss_ppm", "trial", "base_seed_attempt",
        "overhead_cap",
    )
    string_fields = ("phase", "band", "schedule", "loss_seed")
    if (any(type(row[field]) is not int for field in integer_fields) or
            any(not isinstance(row[field], str) for field in string_fields)):
        fail("ledger cell key uses a noncanonical scalar type")
    try:
        K_slot = k_index[row["K"]]
        seed_slot = seed_index[(row["base_seed_attempt"], row["loss_seed"])]
        stratum_slot = stratum_index[(row["schedule"], row["loss_ppm"])]
        width_slot = width_index[row["block_bytes"]]
    except (KeyError, TypeError):
        fail("ledger row is outside the frozen {} domain".format(phase))
    if row["trial"] != seed_slot:
        fail("ledger row swaps a trial index and seed pair")
    expected = {
        key: row[key] for key in contract["recovery"]["cell_key"]
    }
    expected["band"] = _band_for(contract, row["K"])
    if row["phase"] != phase or row["band"] != expected["band"] or \
            row["overhead_cap"] != contract["recovery"]["overhead_cap"]:
        fail("ledger row phase, band, or cap differs from the frozen cell")
    K_count = len(k_index)
    seed_count = len(seed_index)
    stratum_count = len(stratum_index)
    ordinal = (((width_slot * seed_count + seed_slot) * stratum_count +
                stratum_slot) * K_count + K_slot)
    if ordinal >= domain["expected_cells_per_arm"]:
        fail("internal domain ordinal overflow")
    return ordinal, expected


def _parse_canonical_jsonl(
        path: Path, context: str) -> Iterator[Mapping[str, Any]]:
    try:
        source = path.open("rb")
    except OSError as exc:
        fail("cannot open {} {}: {}".format(context, path, exc))
    with source:
        for line_number, line in enumerate(source, 1):
            if len(line) > MAX_JSON_LINE_BYTES:
                fail("{} line {} exceeds {} bytes".format(
                    context, line_number, MAX_JSON_LINE_BYTES))
            if not line.endswith(b"\n") or line in (b"\n", b"\r\n"):
                fail("{} line {} is empty or incomplete".format(
                    context, line_number))
            value = _load_json_bytes(
                line, "{} line {}".format(context, line_number))
            if not isinstance(value, dict):
                fail("{} line {} is not an object".format(context, line_number))
            logical = canonical_json(value).encode("utf-8")
            if line not in (logical + b"\n", logical + b"\r\n"):
                fail("{} line {} is not canonical JSONL".format(
                    context, line_number))
            yield value


def _parse_ledger(path: Path) -> Iterator[Mapping[str, Any]]:
    return _parse_canonical_jsonl(path, "ledger")


def _validate_timing_qualification_controls(
        contract: Mapping[str, Any], phase: str,
        values: Any) -> Tuple[Mapping[str, Any], ...]:
    expected = contract["timing"]["qualification"]["controls"]
    if not isinstance(values, list) or len(values) != len(expected):
        fail("timing qualification must contain exactly two controls")
    result: List[Mapping[str, Any]] = []
    for index, (control_value, policy) in enumerate(zip(values, expected)):
        control = _exact_keys(
            control_value, TIMING_QUALIFICATION_CONTROL_FIELDS,
            "timing qualification control {}".format(index))
        if (control["arm"] != policy["arm"] or
                control["scope"] != policy["scope"]):
            fail("timing qualification substitutes a mandatory control")
        for field in (
                "binary_sha256", "arm_descriptor_sha256",
                "repair_map_sha256"):
            if (not isinstance(control[field], str) or
                    SHA256.fullmatch(control[field]) is None):
                fail("timing qualification control {} is not a SHA-256".
                     format(field))
        if control["arm"] == "wirehair1":
            expected_policy = "not_applicable"
            expect_map = False
        else:
            expected_policy = "raw_base" if phase == "development" else \
                "repair_map"
            expect_map = phase == "final"
        if control["construction_policy"] != expected_policy:
            fail("timing qualification control construction policy changed")
        has_map = control["repair_map_sha256"] != "0" * 64
        if has_map != expect_map:
            fail("timing qualification control repair-map binding changed")
        result.append(control)
    return tuple(result)


def _validate_timing_qualification_outcome(
        row: Mapping[str, Any], arm: str, cap: int) -> str:
    outcome = row[arm + "_outcome"]
    extra = row[arm + "_decoded_extra"]
    if outcome == "success":
        if arm == "wirehair2_head":
            if type(extra) is not int or extra != 4:
                fail("qualified WH2-head solve must succeed exactly at K+4")
        elif type(extra) is not int or not 0 <= extra <= cap:
            fail("qualified Wirehair1 receive has invalid decoded_extra")
    elif outcome == "need_more_at_bound":
        if extra is not None:
            fail("qualification need-more outcome must use decoded_extra=null")
    else:
        fail("only success or need_more_at_bound may enter qualification audit")
    return outcome


def _load_timing_qualification_map(
        contract: Mapping[str, Any], phase: str, path: Path,
        audit_path: Path, expected_sha256: str,
        seal: Any) -> TimingQualification:
    """Load a frozen candidate-blind map and prove every lowest retry."""
    domains = contract["timing"]["domains"]
    domain = domains.get(phase)
    if not isinstance(domain, Mapping):
        fail("unknown timing qualification phase: " + str(phase))
    if (not isinstance(expected_sha256, str) or
            SHA256.fullmatch(expected_sha256) is None):
        fail("expected timing qualification map hash is not a SHA-256")
    value = _load_canonical_json_file(path, "timing qualification map")
    _exact_keys(
        value, TIMING_QUALIFICATION_MAP_FIELDS,
        "timing qualification map")
    actual_map_sha256 = timing_qualification_map_sha256(value)
    if actual_map_sha256 != expected_sha256:
        fail("timing qualification map differs from its frozen SHA-256")
    if (value["schema"] != TIMING_QUALIFICATION_MAP_SCHEMA or
            value["contract_sha256"] != contract_sha256(contract) or
            value["phase"] != phase or
            value["base_domain_sha256"] != domain["base_domain_sha256"] or
            value["entry_kind"] != TIMING_QUALIFICATION_ENTRY_KIND):
        fail("timing qualification map does not bind this contract and phase")
    if (not isinstance(value["source_git_commit"], str) or
            GIT_COMMIT.fullmatch(value["source_git_commit"]) is None):
        fail("timing qualification source commit is not a full Git object ID")
    for field in (
            "qualified_domain_sha256", "qualification_audit_sha256",
            "selected_trace_roster_sha256"):
        if (not isinstance(value[field], str) or
                SHA256.fullmatch(value[field]) is None):
            fail("timing qualification {} is not a SHA-256".format(field))
    controls = _validate_timing_qualification_controls(
        contract, phase, value["controls"])
    retry_offsets_value = value["retry_offsets"]
    expected_cells = domain["expected_cells"]
    if (not isinstance(retry_offsets_value, list) or
            len(retry_offsets_value) != expected_cells):
        fail("timing qualification map needs one retry offset per base cell")
    retry_offsets: List[int] = []
    for ordinal, retry_offset_value in enumerate(retry_offsets_value):
        retry_offset = _integer(
            retry_offset_value,
            "timing qualification retry offset {}".format(ordinal))
        if retry_offset >= TIMING_LOSS_RETRY_ATTEMPTS:
            fail("timing qualification retry offset exceeds uint8 range")
        retry_offsets.append(retry_offset)
    qualified_digest = _timing_domain_sha256_from_offsets(
        contract, phase, retry_offsets)
    if value["qualified_domain_sha256"] != qualified_digest:
        fail("timing qualification qualified-domain SHA-256 is invalid")

    rows = iter(_parse_canonical_jsonl(
        audit_path, "timing qualification audit"))
    audit_digest = _timing_qualification_audit_hasher(contract, phase)
    cap = contract["timing"]["receive_overhead_cap"]
    selected_trace_sha256s: List[str] = []
    for ordinal, (base_cell, selected_offset) in enumerate(zip(
            iter_timing_base_cells(contract, phase), retry_offsets)):
        base_cell_sha256 = sha256_json(base_cell)
        for retry_offset in range(selected_offset + 1):
            try:
                row = next(rows)
            except StopIteration:
                fail("timing qualification audit ends before the selected retry")
            audit_digest.update(canonical_json(row).encode("utf-8"))
            audit_digest.update(b"\n")
            _exact_keys(
                row, TIMING_QUALIFICATION_AUDIT_FIELDS,
                "timing qualification audit row")
            _exact_integer(
                row["ordinal"], ordinal,
                "timing qualification audit ordinal")
            _exact_integer(
                row["loss_retry_offset"], retry_offset,
                "timing qualification audit retry offset")
            if row["base_cell_sha256"] != base_cell_sha256:
                fail("timing qualification audit substitutes its base cell")
            expected_seed = _qualified_timing_loss_seed(
                base_cell["base_loss_seed"], retry_offset)
            if row["loss_seed"] != expected_seed:
                fail("timing qualification audit uses the wrong retry seed")
            if (not isinstance(row["trace_sha256"], str) or
                    SHA256.fullmatch(row["trace_sha256"]) is None):
                fail("timing qualification trace is not a SHA-256")
            outcomes = (
                _validate_timing_qualification_outcome(
                    row, "wirehair2_head", cap),
                _validate_timing_qualification_outcome(
                    row, "wirehair1", cap),
            )
            both_success = outcomes == ("success", "success")
            if retry_offset < selected_offset and both_success:
                fail("timing qualification map does not choose the lowest retry")
            if retry_offset == selected_offset and not both_success:
                fail("timing qualification terminal retry is not successful")
            if retry_offset == selected_offset:
                selected_trace_sha256s.append(row["trace_sha256"])
    try:
        next(rows)
        fail("timing qualification audit contains rows after selected retries")
    except StopIteration:
        pass
    if audit_digest.hexdigest() != value["qualification_audit_sha256"]:
        fail("timing qualification audit differs from its map binding")
    selected_trace_roster_sha256 = timing_selected_trace_roster_sha256(
        selected_trace_sha256s)
    if selected_trace_roster_sha256 != value["selected_trace_roster_sha256"]:
        fail("timing qualification selected trace roster is invalid")
    return seal({
        "contract_sha256": value["contract_sha256"],
        "phase": phase,
        "map_sha256": actual_map_sha256,
        "base_domain_sha256": value["base_domain_sha256"],
        "qualified_domain_sha256": value["qualified_domain_sha256"],
        "source_git_commit": value["source_git_commit"],
        "controls": tuple(
            tuple(control[field]
                  for field in TIMING_QUALIFICATION_CONTROL_ORDER)
            for control in controls),
        "qualification_audit_sha256": value["qualification_audit_sha256"],
        "retry_offsets": tuple(retry_offsets),
        "selected_trace_roster_sha256": selected_trace_roster_sha256,
        "selected_trace_sha256s": tuple(selected_trace_sha256s),
    })


def _make_timing_qualification_loader(seal: Any) -> Any:
    """Capture the registry sealer without publishing it as module authority."""
    def load_timing_qualification_map(
            contract: Mapping[str, Any], phase: str, path: Path,
            audit_path: Path, expected_sha256: str) -> TimingQualification:
        return _load_timing_qualification_map(
            contract, phase, path, audit_path, expected_sha256, seal)
    return load_timing_qualification_map


load_timing_qualification_map = _make_timing_qualification_loader(
    _seal_timing_qualification)
del _seal_timing_qualification


def _trace_manifest_hasher(
        contract: Mapping[str, Any], evidence_kind: str,
        phase: str,
        qualification: Optional[TimingQualification] = None) -> Any:
    domain = contract[evidence_kind]["domains"][phase]
    digest = hashlib.sha256()
    if evidence_kind == "timing":
        if qualification is None:
            fail("timing trace manifests require a validated qualification map")
        validated = _require_timing_qualification(
            contract, phase, qualification)
        digest.update(b"wirehair.wh2.trace-manifest.v2\0")
        digest.update(bytes.fromhex(contract_sha256(contract)))
        digest.update(bytes.fromhex(validated.base_domain_sha256))
        digest.update(bytes.fromhex(validated.qualified_domain_sha256))
        digest.update(bytes.fromhex(validated.map_sha256))
    else:
        digest.update(b"wirehair.wh2.trace-manifest.v1\0")
        digest.update(bytes.fromhex(contract_sha256(contract)))
        digest.update(bytes.fromhex(domain["domain_sha256"]))
    digest.update(phase.encode("utf-8"))
    digest.update(b"\0")
    return digest


def _hash_trace_manifest_row(digest: Any, value: Mapping[str, Any]) -> None:
    digest.update(canonical_json(value).encode("utf-8"))
    digest.update(b"\n")


def trace_manifest_sha256(
        contract: Mapping[str, Any], evidence_kind: str, phase: str,
        path: Path,
        qualification: Optional[TimingQualification] = None) -> str:
    digest = _trace_manifest_hasher(
        contract, evidence_kind, phase, qualification)
    for value in _parse_canonical_jsonl(path, "trace manifest"):
        _hash_trace_manifest_row(digest, value)
    return digest.hexdigest()


def _load_frozen_trace_manifest(
        contract: Mapping[str, Any], evidence_kind: str, phase: str,
        path: Path, expected_sha256: str,
        cells: Iterable[Mapping[str, Any]], count: int,
        context: str,
        qualification: Optional[TimingQualification] = None) -> Sequence[str]:
    digest = _trace_manifest_hasher(
        contract, evidence_kind, phase, qualification)
    traces: List[str] = []
    rows = iter(_parse_canonical_jsonl(path, context))
    for ordinal, cell in enumerate(cells):
        try:
            row = next(rows)
        except StopIteration:
            fail("{} has {} cells, expected {}".format(context, ordinal, count))
        _hash_trace_manifest_row(digest, row)
        if set(row) != TRACE_FIELDS:
            fail("{} row has an unexpected schema".format(context))
        _exact_integer(row["ordinal"], ordinal, context + " ordinal")
        for field in ("cell_sha256", "trace_sha256"):
            if (not isinstance(row[field], str) or
                    SHA256.fullmatch(row[field]) is None):
                fail("{} {} is not a lowercase SHA-256".format(
                    context, field))
        if row["cell_sha256"] != sha256_json(cell):
            fail("{} cell hash does not bind its key".format(context))
        traces.append(row["trace_sha256"])
    if len(traces) != count:
        fail("{} domain yielded {} cells, expected {}".format(
            context, len(traces), count))
    try:
        next(rows)
        fail("{} contains cells beyond the frozen domain".format(context))
    except StopIteration:
        pass
    if digest.hexdigest() != expected_sha256:
        fail("{} differs from the pre-result freeze".format(context))
    return traces


def load_trace_manifest(
        contract: Mapping[str, Any], phase: str, path: Path,
        expected_sha256: str) -> Sequence[str]:
    count = contract["recovery"]["domains"][phase]["expected_cells_per_arm"]
    return _load_frozen_trace_manifest(
        contract, "recovery", phase, path, expected_sha256,
        iter_recovery_cells(contract, phase), count, "trace manifest")


def load_timing_trace_manifest(
        contract: Mapping[str, Any], phase: str, path: Path,
        expected_sha256: str,
        qualification: TimingQualification) -> Sequence[str]:
    count = contract["timing"]["domains"][phase]["expected_cells"]
    traces = _load_frozen_trace_manifest(
        contract, "timing", phase, path, expected_sha256,
        iter_timing_cells(contract, phase, qualification), count,
        "timing trace manifest", qualification)
    if tuple(traces) != qualification.selected_trace_sha256s:
        fail("timing trace manifest does not replay terminal audit traces")
    return traces


def validate_ledger(contract: Mapping[str, Any], phase: str, path: Path,
                    freeze_manifest_path: Path, trace_manifest_path: Path,
                    repair_map_paths: Optional[Mapping[str, Path]] = None,
                    ) -> Dict[str, Any]:
    if repair_map_paths is None:
        repair_map_paths = {}
    freeze = load_freeze_manifest(
        contract, phase, freeze_manifest_path, "recovery")
    arms = tuple(freeze["arm_roster"])
    controls = tuple(contract["selection"]["controls"])
    if freeze["schema"] == RAW_FREEZE_SCHEMA:
        architecture_roles = freeze["architecture_roles"]
        recovery_reference = architecture_roles["recovery_reference"]
        candidate_arms = tuple(architecture_roles["recovery_candidates"])
    else:
        architecture_roles = None
        recovery_reference = None
        candidate_arms = tuple(arm for arm in arms if arm not in controls)
    domain = contract["recovery"]["domains"].get(phase)
    if domain is None:
        fail("unknown recovery phase: " + phase)
    count = domain["expected_cells_per_arm"]
    indexes = _domain_indexes(contract, phase)
    frozen_traces = load_trace_manifest(
        contract, phase, trace_manifest_path, freeze["trace_manifest_sha256"])
    required_map_arms = {
        arm for arm in arms
        if freeze["arms_by_name"][arm]["repair_map_sha256"] != "0" * 64
    }
    if set(repair_map_paths) != required_map_arms:
        fail("repair-map arguments must match the frozen nonzero map roster")
    repair_attempts = {
        arm: load_repair_map(
            contract, repair_map_paths[arm], freeze["arms_by_name"][arm],
            freeze["source_git_commit"],
            freeze["repair_training_trace_manifest_sha256"])
        for arm in required_map_arms
    }
    seen = {arm: bytearray(count) for arm in arms}
    scores = {arm: bytearray(count) for arm in arms}
    status_counts = {arm: defaultdict(int) for arm in arms}
    threshold_counts = {
        arm: {threshold: 0 for threshold in contract["recovery"]["overhead_thresholds"]}
        for arm in arms
    }
    scoped_counts = {arm: defaultdict(lambda: defaultdict(int)) for arm in arms}
    scoped_cells = {arm: defaultdict(int) for arm in arms}
    weak_multiplicity = {arm: defaultdict(int) for arm in arms}
    capped_overhead_sum = {arm: 0 for arm in arms}

    for row in _parse_ledger(path):
        if not isinstance(row, dict):
            fail("ledger row has an unexpected schema")
        arm = row.get("arm")
        if not isinstance(arm, str):
            fail("ledger arm name must be a string")
        if arm not in seen:
            fail("ledger contains an arm outside the frozen roster")
        frozen_arm = freeze["arms_by_name"][arm]
        raw_arm = freeze["schema"] in (
            LEGACY_RAW_FREEZE_SCHEMA, RAW_FREEZE_SCHEMA) and \
            frozen_arm["construction_seed_basis"] == \
            RAW_CONSTRUCTION_SEED_BASIS
        if raw_arm:
            expected_fields = RAW_RECOVERY_RECORD_FIELDS if \
                freeze["schema"] == RAW_FREEZE_SCHEMA else \
                LEGACY_RAW_RECOVERY_RECORD_FIELDS
        else:
            expected_fields = LEDGER_FIELDS
        if set(row) != expected_fields:
            fail("ledger row has an unexpected schema")
        for field in ("cell_sha256", "trace_sha256", "binary_sha256",
                      "arm_descriptor_sha256", "realized_construction_sha256",
                      "repair_map_sha256"):
            if not isinstance(row[field], str) or SHA256.fullmatch(row[field]) is None:
                fail("ledger {} is not a lowercase SHA-256".format(field))
        ordinal, cell = _cell_ordinal(contract, phase, row, indexes)
        if row["cell_sha256"] != sha256_json(cell):
            fail("ledger cell hash does not bind its arm-free key")
        if seen[arm][ordinal]:
            fail("duplicate ledger cell for arm {} at ordinal {}".format(arm, ordinal))
        seen[arm][ordinal] = 1
        if row["trace_sha256"] != frozen_traces[ordinal]:
            fail("ledger trace differs from the pre-result trace manifest")
        for field in ("binary_sha256", "arm_descriptor_sha256",
                      "repair_map_sha256"):
            if row[field] != frozen_arm[field]:
                fail("ledger {} differs from frozen arm {}".format(field, arm))

        attempt = row["construction_attempt"]
        if (type(attempt) is not int or
                not 0 <= attempt < contract["recovery"]["max_construction_attempts"]):
            fail("ledger construction attempt is outside [0,255]")
        expected_attempt = expected_construction_attempt(
            frozen_arm, row["K"], row["base_seed_attempt"], repair_attempts)
        if attempt != expected_attempt:
            fail("ledger construction attempt differs from the frozen seed policy")
        if raw_arm:
            construction_fields = RAW_CONSTRUCTION_FIELDS if \
                freeze["schema"] == RAW_FREEZE_SCHEMA else \
                LEGACY_RAW_CONSTRUCTION_FIELDS
            raw_fields = {
                field: row[field] for field in construction_fields
            }
            _validate_raw_construction_fields(
                raw_fields, "raw ledger row",
                freeze["schema"] == RAW_FREEZE_SCHEMA)
            _validate_raw_descriptor_structure(
                arm, frozen_arm["arm_descriptor_sha256"], row["K"], raw_fields,
                "raw ledger row")
            if (row["construction_seed_basis"] !=
                    frozen_arm["construction_seed_basis"] or
                    row["seed_schedule_sha256"] !=
                    frozen_arm["seed_schedule_sha256"]):
                fail("raw ledger seed policy differs from the frozen arm")
            if (row["precode_attempt"] != attempt or
                    row["packet_attempt"] != attempt):
                fail("raw ledger attempts differ from the frozen paired attempt")
            if freeze["schema"] == RAW_FREEZE_SCHEMA:
                expected_realized = raw_realized_construction_sha256(
                    frozen_arm["codec"], arm,
                    frozen_arm["arm_descriptor_sha256"], row["K"],
                    row["block_bytes"], raw_fields)
            else:
                expected_realized = legacy_raw_realized_construction_sha256(
                    frozen_arm["codec"], arm,
                    frozen_arm["arm_descriptor_sha256"], row["K"],
                    row["block_bytes"], raw_fields)
        else:
            expected_realized = generic_realized_construction_sha256(
                frozen_arm["codec"], frozen_arm["arm_descriptor_sha256"],
                row["K"], row["block_bytes"], attempt)
        if row["realized_construction_sha256"] != expected_realized:
            fail("ledger realized construction does not match its frozen descriptor")

        outcome = row["outcome"]
        extra = row["decoded_extra"]
        cap = contract["recovery"]["overhead_cap"]
        if outcome not in SCORED_OUTCOMES:
            fail("fatal, internal, or unknown outcomes abort rather than score")
        if outcome == "success":
            if type(extra) is not int or not 0 <= extra <= cap:
                fail("successful ledger row has invalid decoded_extra")
            score = extra
        else:
            if extra is not None:
                fail("failed ledger row must use decoded_extra=null")
            score = cap + 1
        scores[arm][ordinal] = score
        status_counts[arm][outcome] += 1
        capped_overhead_sum[arm] += min(score, cap)
        scope = (
            row["band"], row["block_bytes"], row["schedule"], row["loss_ppm"])
        scoped_cells[arm][scope] += 1
        for threshold in threshold_counts[arm]:
            failed = score > threshold
            threshold_counts[arm][threshold] += failed
            scoped_counts[arm][scope][threshold] += failed
        if score > 0:
            weak_multiplicity[arm][(
                row["K"], row["base_seed_attempt"], row["block_bytes"])] += 1

    for arm in arms:
        present = sum(seen[arm])
        if present != count:
            fail("arm {} has {} cells, expected {} (no exclusions allowed)".format(
                arm, present, count))
    summaries = {}
    target_ppm = contract["goal"]["primary_failure_threshold_ppm"]
    for arm in arms:
        primary_failures = threshold_counts[arm][0]
        scope_summaries = {}
        for scope, failures in sorted(scoped_counts[arm].items()):
            scope_cells = scoped_cells[arm][scope]
            primary = failures[0]
            scope_summaries["{}|{}|{}|{}".format(*scope)] = {
                "cells": scope_cells,
                "failure_by_overhead": {
                    str(key): value for key, value in failures.items()
                },
                "one_percent_overhead0_pass":
                    primary * 1000000 <= target_ppm * scope_cells,
            }
        summaries[arm] = {
            "cells": count,
            "failure_by_overhead": {
                str(key): value for key, value in threshold_counts[arm].items()
            },
            "failure_ppm_overhead0": (primary_failures * 1000000) // count,
            "one_percent_overall_pass": primary_failures * 1000000 <= target_ppm * count,
            "status_counts": dict(sorted(status_counts[arm].items())),
            "capped_overhead_sum": capped_overhead_sum[arm],
            "weak_seed_units": len(weak_multiplicity[arm]),
            "weak_seed_multiplicity_histogram": _histogram(weak_multiplicity[arm].values()),
            "weak_seed_unit_details": [
                {
                    "K": key[0], "base_seed_attempt": key[1],
                    "block_bytes": key[2],
                    "overhead0_failure_multiplicity": multiplicity,
                }
                for key, multiplicity in sorted(weak_multiplicity[arm].items())
            ],
            "all_band_strata_one_percent_pass": all(
                value["one_percent_overhead0_pass"]
                for value in scope_summaries.values()
            ),
            "band_stratum": scope_summaries,
        }
        construct_failed = status_counts[arm].get("construct_failed", 0)
        unsupported = status_counts[arm].get("unsupported", 0)
        summaries[arm]["construct_failed"] = construct_failed
        summaries[arm]["unsupported"] = unsupported
        if phase in ("final_repaired", "final_validation"):
            summaries[arm]["phase_recovery_gate_pass"] = (
                len(weak_multiplicity[arm]) == 0 and
                construct_failed == 0 and unsupported == 0)
        elif phase == "final_raw":
            summaries[arm]["phase_recovery_gate_pass"] = (
                summaries[arm]["one_percent_overall_pass"] and
                summaries[arm]["all_band_strata_one_percent_pass"] and
                unsupported == 0)
        elif phase == "cross_width_validation":
            summaries[arm]["phase_recovery_gate_pass"] = (
                summaries[arm]["one_percent_overall_pass"] and
                summaries[arm]["all_band_strata_one_percent_pass"] and
                construct_failed == 0 and unsupported == 0)
        else:
            summaries[arm]["phase_recovery_gate_pass"] = unsupported == 0

    comparisons = {}
    margin_ppm = contract["selection"]["noninferiority_margin_ppm"]
    mandatory_controls_supported = all(
        summaries[control]["unsupported"] == 0 for control in controls)
    for arm in candidate_arms:
        control_comparisons = {}
        comparison_arms = controls if recovery_reference is None else \
            (recovery_reference,) + controls
        for control in comparison_arms:
            repairs = introductions = shared_failures = 0
            for left, right in zip(scores[control], scores[arm]):
                left_failed = left > 0
                right_failed = right > 0
                repairs += left_failed and not right_failed
                introductions += not left_failed and right_failed
                shared_failures += left_failed and right_failed
            allowed = max(
                1, (margin_ppm * count + 1000000 - 1) // 1000000)
            overall_tail = {
                str(threshold): {
                    "control_failures": threshold_counts[control][threshold],
                    "candidate_failures": threshold_counts[arm][threshold],
                    "allowed_excess_failures": allowed,
                    "noninferior": threshold_counts[arm][threshold] <=
                        threshold_counts[control][threshold] + allowed,
                }
                for threshold in contract["recovery"]["overhead_thresholds"]
            }
            scope_tail = {}
            for scope in sorted(scoped_cells[arm]):
                scope_count = scoped_cells[arm][scope]
                scope_allowed = max(
                    1, (margin_ppm * scope_count + 1000000 - 1) // 1000000)
                scope_tail["{}|{}|{}|{}".format(*scope)] = {
                    str(threshold): {
                        "control_failures": scoped_counts[control][scope][threshold],
                        "candidate_failures": scoped_counts[arm][scope][threshold],
                        "allowed_excess_failures": scope_allowed,
                        "noninferior": scoped_counts[arm][scope][threshold] <=
                            scoped_counts[control][scope][threshold] + scope_allowed,
                    }
                    for threshold in contract["recovery"]["overhead_thresholds"]
                }
            control_comparisons[control] = {
                "overhead0_repairs": repairs,
                "overhead0_introductions": introductions,
                "overhead0_shared_failures": shared_failures,
                "overall_tail": overall_tail,
                "band_stratum_tail": scope_tail,
                "all_noninferiority_gates_pass": (
                    all(value["noninferior"] for value in overall_tail.values()) and
                    all(value["noninferior"] for tail in scope_tail.values()
                        for value in tail.values())
                ),
            }
        if recovery_reference is None:
            all_controls_noninferior = all(
                value["all_noninferiority_gates_pass"]
                for value in control_comparisons.values()) and \
                summaries[arm]["unsupported"] == 0 and \
                mandatory_controls_supported
            comparisons[arm] = {
                "controls": control_comparisons,
                "all_controls_noninferior": all_controls_noninferior,
                "architecture_eligible": (
                    all_controls_noninferior and
                    summaries[arm]["phase_recovery_gate_pass"]),
            }
        else:
            reference_comparison = control_comparisons[recovery_reference]
            reference_noninferior = \
                reference_comparison["all_noninferiority_gates_pass"] and \
                summaries[arm]["unsupported"] == 0 and \
                summaries[recovery_reference]["unsupported"] == 0
            comparisons[arm] = {
                "recovery_reference": recovery_reference,
                "reference_comparison": reference_comparison,
                "descriptive_controls": {
                    control: control_comparisons[control]
                    for control in controls
                },
                "reference_noninferior": reference_noninferior,
                "architecture_eligible": (
                    reference_noninferior and
                    summaries[arm]["phase_recovery_gate_pass"] and
                    summaries[recovery_reference][
                        "phase_recovery_gate_pass"]),
            }
    result = {
        "schema": SCHEMA + ".ledger-summary.v1",
        "phase": phase,
        "contract_sha256": contract_sha256(contract),
        "freeze_manifest_sha256": freeze_manifest_sha256(freeze),
        "architecture_artifact_sha256":
            architecture_artifact_sha256(freeze),
        "arm_artifacts": frozen_arm_artifacts(freeze),
        "domain_sha256": domain["domain_sha256"],
        "excluded_cells": 0,
        "mandatory_controls_supported": mandatory_controls_supported,
        "arms": summaries,
        "comparisons": comparisons,
    }
    if freeze["schema"] == LEGACY_RAW_FREEZE_SCHEMA:
        result["schema"] = SCHEMA + ".ledger-summary.v2"
        result["freeze_schema"] = LEGACY_RAW_FREEZE_SCHEMA
    elif freeze["schema"] == RAW_FREEZE_SCHEMA:
        result["schema"] = SCHEMA + ".ledger-summary.v3"
        result["freeze_schema"] = RAW_FREEZE_SCHEMA
        result["architecture_roles"] = architecture_roles
        result["timing_proxy_witness_sha256"] = \
            freeze["timing_proxy_witness_sha256"]
        result["work_rank_summary_sha256"] = \
            freeze["work_rank_summary_sha256"]
        result["work_rank_result_stream_sha256"] = \
            freeze["work_rank_result_stream_sha256"]
        result["work_rank_domain_sha256"] = \
            freeze["work_rank_domain_sha256"]
    return result


def _timing_outcome(
        row: Mapping[str, Any], side: str, scope: str) -> Tuple[str, Any]:
    outcome = row[side + "_outcome"]
    extra = row[side + "_decoded_extra"]
    if outcome not in SCORED_OUTCOMES:
        fail("fatal, internal, or unknown timing outcomes abort")
    if outcome == "success":
        if scope == "encoder_init_plus_first_K_symbols":
            if extra is not None:
                fail("successful encoder timing must use decoded_extra=null")
        elif scope == "decoder_solve":
            if type(extra) is not int or extra != 4:
                fail("isolated solve timing must receipt the fixed K+4 prefix")
        elif (type(extra) is not int or
              not 0 <= extra <= row["receive_overhead_cap"]):
            fail("successful decoder timing has invalid decoded_extra")
    elif extra is not None:
        fail("failed timing side must use decoded_extra=null")
    return outcome, extra


def _timing_identity(
        row: Mapping[str, Any], side: str, arm: str,
        frozen_arm: Mapping[str, Any], repair_attempts: Mapping[str, Mapping[int, int]],
        ) -> Tuple[Any, ...]:
    hash_fields = (
        "binary_sha256", "arm_descriptor_sha256",
        "realized_construction_sha256", "repair_map_sha256",
    )
    for field in hash_fields:
        value = row[side + "_" + field]
        if not isinstance(value, str) or SHA256.fullmatch(value) is None:
            fail("timing {}_{} is not a lowercase SHA-256".format(side, field))
    for field in ("binary_sha256", "arm_descriptor_sha256",
                  "repair_map_sha256"):
        if row[side + "_" + field] != frozen_arm[field]:
            fail("timing {} {} differs from the frozen arm".format(side, field))
    attempt = row[side + "_construction_attempt"]
    if type(attempt) is not int:
        fail("timing construction attempt is not an integer")
    expected_attempt = expected_construction_attempt(
        frozen_arm, row["K"], row["base_seed_attempt"], repair_attempts)
    if attempt != expected_attempt:
        fail("timing construction attempt differs from the frozen seed policy")
    expected_realized = generic_realized_construction_sha256(
        frozen_arm["codec"], frozen_arm["arm_descriptor_sha256"],
        row["K"], row["block_bytes"], attempt)
    if row[side + "_realized_construction_sha256"] != expected_realized:
        fail("timing realized construction differs from its frozen descriptor")
    return (
        arm, frozen_arm["binary_sha256"], frozen_arm["arm_descriptor_sha256"],
        attempt, expected_realized, frozen_arm["repair_map_sha256"],
    )


def _timing_confidence_interval(
        contract: Mapping[str, Any], values: Sequence[float],
        expected_count: int) -> Mapping[str, float]:
    if len(values) != expected_count or expected_count < 2:
        fail("timing statistics lack the complete replicate domain")
    mean = math.fsum(values) / expected_count
    variance = math.fsum(
        (value - mean) * (value - mean) for value in values) / \
        (expected_count - 1)
    critical = contract["timing"]["statistics"][
        "t_critical_by_independent_rounds"].get(str(expected_count))
    if type(critical) is not float:
        fail("timing contract lacks the exact Student-t critical value")
    half_width = critical * math.sqrt(variance / expected_count)
    return {
        "mean_log_ratio": mean,
        "lower_95": mean - half_width,
        "upper_95": mean + half_width,
    }


def _timing_group_decision(
        contract: Mapping[str, Any], ab: Mapping[str, Any],
        left_aa: Mapping[str, Any],
        right_aa: Mapping[str, Any]) -> Mapping[str, Any]:
    if not (ab.get("selectable") and left_aa.get("selectable") and
            right_aa.get("selectable")):
        return {
            "selectable": False, "status": "non_selectable",
            "reason": "failed_or_incomplete_panel",
        }
    floor = math.log1p(
        contract["timing"]["practical_margin_ppm"] / 1000000.0)
    if any(not (-floor < aa["lower_95"] and aa["upper_95"] < floor)
           for aa in (left_aa, right_aa)):
        return {
            "selectable": False, "status": "non_selectable",
            "reason": "aa_repeatability_threshold",
            "effective_floor": floor,
        }
    if ab["upper_95"] < -floor:
        status = "faster"
    elif ab["lower_95"] > floor:
        status = "resolved_slower"
    elif ab["upper_95"] < floor:
        status = "noninferior"
    else:
        status = "unresolved"
    return {
        "selectable": True,
        "status": status,
        "effective_floor": floor,
        **ab,
    }


def _timing_speed_gate(
        phase: str, candidate_decisions: Mapping[str, Any]) -> bool:
    if phase == "development":
        all_groups = [
            decision for panel in candidate_decisions.values()
            for decision in panel["by_band_width"].values()
        ] + [
            panel["aggregate"] for panel in candidate_decisions.values()
        ]
        return all(
            value["status"] not in ("resolved_slower", "non_selectable")
            for value in all_groups)
    if phase != "final":
        fail("unknown timing speed-gate phase: " + phase)
    solve = candidate_decisions["decoder_solve|wirehair2_head"]
    receive = candidate_decisions["receive_to_success|wirehair1"]
    encoder_head = candidate_decisions[
        "encoder_init_plus_first_K_symbols|wirehair2_head"]
    encoder_wh1 = candidate_decisions[
        "encoder_init_plus_first_K_symbols|wirehair1"]
    return (
        all(value["status"] == "faster"
            for value in solve["by_band_width"].values()) and
        all(value["status"] == "faster"
            for value in receive["by_band_width"].values()) and
        all(value["status"] in ("faster", "noninferior")
            for panel in (encoder_head, encoder_wh1)
            for value in panel["by_band_width"].values()) and
        encoder_head["aggregate"]["status"] == "faster" and
        encoder_wh1["aggregate"]["status"] == "faster")


def validate_timing_receipt(
        contract: Mapping[str, Any], phase: str, path: Path,
        freeze_manifest_path: Path, trace_manifest_path: Path,
        timing_qualification_map_path: Path,
        timing_qualification_audit_path: Path,
        timing_qualification_map_digest: str,
        repair_map_paths: Optional[Mapping[str, Path]] = None,
        ) -> Dict[str, Any]:
    if repair_map_paths is None:
        repair_map_paths = {}
    qualification = load_timing_qualification_map(
        contract, phase, timing_qualification_map_path,
        timing_qualification_audit_path, timing_qualification_map_digest)
    freeze = load_freeze_manifest(
        contract, phase, freeze_manifest_path, "timing", qualification)
    arms = tuple(freeze["arm_roster"])
    controls = tuple(contract["selection"]["controls"])
    candidates = tuple(arm for arm in arms if arm not in controls)
    domain = contract["timing"]["domains"].get(phase)
    if domain is None:
        fail("unknown timing phase: " + phase)
    cell_count = domain["expected_cells"]
    cells = _timing_cell_indexes(contract, phase, qualification)
    traces = load_timing_trace_manifest(
        contract, phase, trace_manifest_path, freeze["trace_manifest_sha256"],
        qualification)
    required_map_arms = {
        arm for arm in arms
        if freeze["arms_by_name"][arm]["repair_map_sha256"] != "0" * 64
    }
    if set(repair_map_paths) != required_map_arms:
        fail("timing repair-map arguments must match the frozen map roster")
    repair_attempts = {
        arm: load_repair_map(
            contract, repair_map_paths[arm], freeze["arms_by_name"][arm],
            freeze["source_git_commit"],
            freeze["repair_training_trace_manifest_sha256"])
        for arm in required_map_arms
    }
    panels = timing_panels(contract, arms)
    panel_indexes = {
        canonical_json(panel): index for index, panel in enumerate(panels)
    }
    seen = bytearray(cell_count * len(panels))
    measurements: Dict[str, Any] = {
        canonical_json(panel): defaultdict(list) for panel in panels
    }
    failed_groups: Dict[str, Any] = {
        canonical_json(panel): set() for panel in panels
    }
    failed_aggregate = {canonical_json(panel): False for panel in panels}
    construction_by_arm_cell: Dict[Tuple[str, int], Tuple[Any, ...]] = {}
    outcome_by_arm_cell_scope: Dict[
        Tuple[str, int, str], Tuple[str, Any]] = {}

    for row in _parse_canonical_jsonl(path, "timing receipt"):
        if set(row) != TIMING_RECEIPT_FIELDS:
            fail("timing receipt row has an unexpected schema")
        for field in ("cell_sha256", "trace_sha256"):
            if not isinstance(row[field], str) or SHA256.fullmatch(row[field]) is None:
                fail("timing receipt {} is not a lowercase SHA-256".format(field))
        cell_ordinal, cell = _timing_cell_ordinal(contract, phase, row, cells)
        if row["cell_sha256"] != sha256_json(cell):
            fail("timing receipt cell hash does not bind its key")
        if row["trace_sha256"] != traces[cell_ordinal]:
            fail("timing receipt trace differs from the pre-result manifest")
        panel = {
            "panel_kind": row["panel_kind"], "scope": row["scope"],
            "left_arm": row["left_arm"], "right_arm": row["right_arm"],
        }
        panel_key = canonical_json(panel)
        if panel_key not in panel_indexes:
            fail("timing receipt has an undeclared or side-swapped panel")
        receipt_ordinal = cell_ordinal * len(panels) + panel_indexes[panel_key]
        if seen[receipt_ordinal]:
            fail("duplicate timing receipt cell/panel")
        seen[receipt_ordinal] = 1
        expected_order = timing_order(panel, row["replicate"])
        if row["order"] != expected_order:
            fail("timing receipt order differs from frozen counterbalancing")

        side_outcomes = []
        side_identities = []
        for side in ("left", "right"):
            arm = row[side + "_arm"]
            if arm not in freeze["arms_by_name"]:
                fail("timing receipt side is outside the frozen roster")
            frozen_arm = freeze["arms_by_name"][arm]
            identity = _timing_identity(
                row, side, arm, frozen_arm, repair_attempts)
            outcome = _timing_outcome(row, side, row["scope"])
            identity_key = (arm, cell_ordinal)
            prior_identity = construction_by_arm_cell.setdefault(
                identity_key, identity)
            if prior_identity != identity:
                fail("timing arm construction changes between panels")
            outcome_key = (arm, cell_ordinal, row["scope"])
            prior_outcome = outcome_by_arm_cell_scope.setdefault(
                outcome_key, outcome)
            if prior_outcome != outcome:
                fail("timing arm outcome changes between equivalent panels")
            side_identities.append(identity)
            side_outcomes.append(outcome)
        if row["panel_kind"] == "AA" and (
                side_identities[0] != side_identities[1] or
                side_outcomes[0] != side_outcomes[1]):
            fail("A/A timing sides must have identical identity and outcome")

        elapsed = row["elapsed_ns"]
        both_success = all(outcome[0] == "success" for outcome in side_outcomes)
        if not isinstance(elapsed, list) or len(elapsed) != 8:
            fail("timing elapsed_ns must contain exactly eight slots")
        if both_success:
            if any(type(value) is not int or not 1 <= value < (1 << 63)
                   for value in elapsed):
                fail("successful timing slots must be positive int63 nanoseconds")
            logs = [math.log(value) for value in elapsed]
            opposite_order = "BAAB" if expected_order == "ABBA" else "ABBA"
            primary_contrast = _timing_four_slot_log_contrast(
                logs, 0, expected_order)
            opposite_contrast = _timing_four_slot_log_contrast(
                logs, 4, opposite_order)
            log_ratio = (primary_contrast + opposite_contrast) / 2.0
            block_phase_difference = primary_contrast - opposite_contrast
            group = (row["band"], row["block_bytes"], row["replicate"])
            measurements[panel_key][group].append(
                (log_ratio, block_phase_difference))
        else:
            if any(value is not None for value in elapsed):
                fail("failed timing panels must use eight null elapsed slots")
            group_scope = (row["band"], row["block_bytes"])
            failed_groups[panel_key].add(group_scope)
            failed_aggregate[panel_key] = True

    expected_rows = cell_count * len(panels)
    if sum(seen) != expected_rows:
        fail("timing receipt has {} rows, expected {}".format(
            sum(seen), expected_rows))

    repetitions = domain["paired_repetitions"]
    independent_rounds = domain["independent_rounds"]
    lanes_per_round = contract["timing"]["execution_geometry"][
        "jobs_per_wave"]
    if repetitions != independent_rounds * lanes_per_round:
        fail("timing repetitions do not form complete independent rounds")
    K_values = list(_k_values(contract, domain["k_set"]))
    K_count_by_band: Dict[str, int] = defaultdict(int)
    for K in K_values:
        K_count_by_band[_band_for(contract, K)] += 1
    panel_statistics: Dict[str, Any] = {}
    for panel in panels:
        panel_key = canonical_json(panel)
        by_band_width = {}
        for band in (value[0] for value in EXPECTED_BANDS):
            expected_K_count = K_count_by_band[band]
            for width in domain["block_bytes"]:
                label = "{}|{}".format(band, width)
                if (band, width) in failed_groups[panel_key]:
                    by_band_width[label] = {"selectable": False}
                    continue
                round_means = []
                block_phase_round_means = []
                for independent_round in range(independent_rounds):
                    lane_means = []
                    block_phase_lane_means = []
                    for lane in range(lanes_per_round):
                        replicate = independent_round * lanes_per_round + lane
                        values = measurements[panel_key][(
                            band, width, replicate)]
                        if len(values) != expected_K_count:
                            fail("timing band/width lane has an incomplete K "
                                 "cohort")
                        lane_means.append(math.fsum(
                            value[0] for value in values) / len(values))
                        block_phase_lane_means.append(math.fsum(
                            value[1] for value in values) / len(values))
                    round_means.append(
                        math.fsum(lane_means) / lanes_per_round)
                    block_phase_round_means.append(
                        math.fsum(block_phase_lane_means) / lanes_per_round)
                by_band_width[label] = {
                    "selectable": True,
                    **_timing_confidence_interval(
                        contract, round_means, independent_rounds),
                    "block_phase_difference": _timing_confidence_interval(
                        contract, block_phase_round_means,
                        independent_rounds),
                }
        if failed_aggregate[panel_key]:
            aggregate = {"selectable": False}
        else:
            aggregate_round_means = []
            block_phase_round_means = []
            expected_per_lane = len(K_values) * len(domain["block_bytes"])
            for independent_round in range(independent_rounds):
                lane_means = []
                block_phase_lane_means = []
                for lane in range(lanes_per_round):
                    replicate = independent_round * lanes_per_round + lane
                    values = []
                    for band in (value[0] for value in EXPECTED_BANDS):
                        for width in domain["block_bytes"]:
                            values.extend(measurements[panel_key][(
                                band, width, replicate)])
                    if len(values) != expected_per_lane:
                        fail("aggregate timing lane has an incomplete cell "
                             "cohort")
                    lane_means.append(math.fsum(
                        value[0] for value in values) / len(values))
                    block_phase_lane_means.append(math.fsum(
                        value[1] for value in values) / len(values))
                aggregate_round_means.append(
                    math.fsum(lane_means) / lanes_per_round)
                block_phase_round_means.append(
                    math.fsum(block_phase_lane_means) / lanes_per_round)
            aggregate = {
                "selectable": True,
                **_timing_confidence_interval(
                    contract, aggregate_round_means, independent_rounds),
                "block_phase_difference": _timing_confidence_interval(
                    contract, block_phase_round_means,
                    independent_rounds),
            }
        panel_statistics[panel_key] = {
            "panel": panel,
            "by_band_width": by_band_width,
            "aggregate": aggregate,
        }

    decisions: Dict[str, Any] = {}
    for candidate in candidates:
        candidate_decisions = {}
        for value in contract["timing"]["panel_protocol"]["candidate_ab"]:
            panel = {
                "panel_kind": "AB", "scope": value["scope"],
                "left_arm": candidate, "right_arm": value["control"],
            }
            left_aa = {
                "panel_kind": "AA", "scope": value["scope"],
                "left_arm": candidate, "right_arm": candidate,
            }
            right_aa = {
                "panel_kind": "AA", "scope": value["scope"],
                "left_arm": value["control"], "right_arm": value["control"],
            }
            ab_stats = panel_statistics[canonical_json(panel)]
            left_stats = panel_statistics[canonical_json(left_aa)]
            right_stats = panel_statistics[canonical_json(right_aa)]
            group_decisions = {
                label: _timing_group_decision(
                    contract, stats, left_stats["by_band_width"][label],
                    right_stats["by_band_width"][label])
                for label, stats in ab_stats["by_band_width"].items()
            }
            aggregate_decision = _timing_group_decision(
                contract, ab_stats["aggregate"], left_stats["aggregate"],
                right_stats["aggregate"])
            key = "{}|{}".format(value["scope"], value["control"])
            candidate_decisions[key] = {
                "by_band_width": group_decisions,
                "aggregate": aggregate_decision,
            }
        decisions[candidate] = {
            "panels": candidate_decisions,
            "phase_speed_gate_pass": _timing_speed_gate(
                phase, candidate_decisions),
        }

    return {
        "schema": SCHEMA + ".timing-summary.v1",
        "phase": phase,
        "contract_sha256": contract_sha256(contract),
        "base_domain_sha256": domain["base_domain_sha256"],
        "domain_sha256": qualification.qualified_domain_sha256,
        "timing_qualification_map_sha256": qualification.map_sha256,
        "freeze_manifest_sha256": freeze_manifest_sha256(freeze),
        "architecture_artifact_sha256":
            architecture_artifact_sha256(freeze),
        "arm_artifacts": frozen_arm_artifacts(freeze),
        "rows": expected_rows,
        "excluded_cells": 0,
        "panel_statistics": panel_statistics,
        "candidates": decisions,
    }


def select_development_architecture(
        contract: Mapping[str, Any], recovery: Mapping[str, Any],
        timing: Mapping[str, Any],
        timing_qualification: TimingQualification) -> Mapping[str, Any]:
    """Apply the frozen architecture ordering to validated development summaries."""
    contract_digest = contract_sha256(contract)
    recovery_domain = contract["recovery"]["domains"]["development"]
    timing_domain = contract["timing"]["domains"]["development"]
    qualification = _require_timing_qualification(
        contract, "development", timing_qualification)
    if (recovery.get("schema") != SCHEMA + ".ledger-summary.v2" or
            timing.get("schema") != SCHEMA + ".timing-summary.v1" or
            recovery.get("phase") != "development" or
            recovery.get("freeze_schema") != RAW_FREEZE_SCHEMA or
            timing.get("phase") != "development" or
            recovery.get("contract_sha256") != contract_digest or
            timing.get("contract_sha256") != contract_digest or
            recovery.get("domain_sha256") != recovery_domain["domain_sha256"] or
            timing.get("base_domain_sha256") !=
                timing_domain["base_domain_sha256"] or
            timing.get("domain_sha256") !=
                qualification.qualified_domain_sha256 or
            timing.get("timing_qualification_map_sha256") !=
                qualification.map_sha256 or
            type(recovery.get("excluded_cells")) is not int or
            recovery.get("excluded_cells") != 0 or
            type(timing.get("excluded_cells")) is not int or
            timing.get("excluded_cells") != 0):
        fail("architecture selection requires development validator summaries")
    recovery_freeze = recovery.get("freeze_manifest_sha256")
    timing_freeze = timing.get("freeze_manifest_sha256")
    if (not isinstance(recovery_freeze, str) or
            SHA256.fullmatch(recovery_freeze) is None or
            not isinstance(timing_freeze, str) or
            SHA256.fullmatch(timing_freeze) is None):
        fail("development summaries lack frozen manifest identities")
    recovery_artifacts = recovery.get("architecture_artifact_sha256")
    if (not isinstance(recovery_artifacts, str) or
            SHA256.fullmatch(recovery_artifacts) is None or
            recovery_artifacts != timing.get("architecture_artifact_sha256")):
        fail("development recovery and timing substitute architecture artifacts")
    arm_artifacts = recovery.get("arm_artifacts")
    if (not isinstance(arm_artifacts, dict) or
            arm_artifacts != timing.get("arm_artifacts")):
        fail("development summaries substitute per-arm artifacts")
    recovery_comparisons = recovery.get("comparisons")
    timing_candidates = timing.get("candidates")
    recovery_arms = recovery.get("arms")
    if (not isinstance(recovery_comparisons, dict) or
            not isinstance(timing_candidates, dict) or
            not isinstance(recovery_arms, dict) or
            set(recovery_comparisons) != set(timing_candidates)):
        fail("development summaries disagree on the candidate roster")
    candidates = sorted(recovery_comparisons)
    if not candidates:
        fail("development summaries contain no candidate")
    controls = tuple(contract["selection"]["controls"])
    if set(candidates) & set(controls):
        fail("reserved controls cannot enter the candidate roster")
    if set(recovery_arms) != set(controls) | set(candidates):
        fail("development recovery summary has the wrong arm roster")
    if set(arm_artifacts) != set(recovery_arms):
        fail("development per-arm artifacts have the wrong roster")
    for arm, artifact_value in arm_artifacts.items():
        artifact = _exact_keys(artifact_value, (
            "codec", "binary_sha256", "arm_descriptor_sha256",
        ), "development arm artifact " + arm)
        if (not isinstance(artifact["codec"], str) or
                any(not isinstance(artifact[field], str) or
                    SHA256.fullmatch(artifact[field]) is None
                    for field in ("binary_sha256", "arm_descriptor_sha256"))):
            fail("development per-arm artifact is malformed")
    panel_count = len(contract["timing"]["panel_protocol"]["control_aa"]) + \
        len(candidates) * (
            len(contract["timing"]["panel_protocol"]["candidate_aa_scopes"]) +
            len(contract["timing"]["panel_protocol"]["candidate_ab"]))
    expected_timing_rows = timing_domain["expected_cells"] * panel_count
    if (type(timing.get("rows")) is not int or
            timing.get("rows") != expected_timing_rows or
            type(recovery.get("mandatory_controls_supported")) is not bool):
        fail("development summaries have the wrong frozen cardinality")

    eligible = []
    failures: Dict[str, Mapping[str, int]] = {}
    thresholds = contract["recovery"]["overhead_thresholds"]
    for arm in sorted(recovery_arms):
        arm_summary = recovery_arms.get(arm)
        if not isinstance(arm_summary, dict):
            fail("development arm summary is malformed")
        arm_cells = arm_summary.get("cells")
        tail = arm_summary.get("failure_by_overhead")
        if (type(arm_cells) is not int or
                arm_cells != recovery_domain["expected_cells_per_arm"] or
                not isinstance(tail, dict) or set(tail) != {
                    str(threshold) for threshold in thresholds} or
                any(type(tail.get(str(threshold))) is not int or
                    not 0 <= tail[str(threshold)] <= arm_cells
                    for threshold in thresholds)):
            fail("development recovery totals are malformed")
        ordered_tail = [tail[str(threshold)] for threshold in thresholds]
        if any(left < right for left, right in zip(
                ordered_tail, ordered_tail[1:])):
            fail("development recovery tails are not nested")
        failures[arm] = tail

    for arm in candidates:
        comparison = recovery_comparisons[arm]
        timing_summary = timing_candidates[arm]
        if (not isinstance(comparison, dict) or
                not isinstance(timing_summary, dict) or
                type(comparison.get("architecture_eligible")) is not bool or
                type(timing_summary.get("phase_speed_gate_pass")) is not bool or
                not isinstance(comparison.get("controls"), dict) or
                set(comparison["controls"]) != set(controls)):
            fail("development candidate summary is malformed")
        if (comparison.get("architecture_eligible") is True and
                timing_summary.get("phase_speed_gate_pass") is True and
                recovery["mandatory_controls_supported"] is True):
            eligible.append(arm)

    bindings = {
        "contract_sha256": contract_digest,
        "recovery_domain_sha256": recovery_domain["domain_sha256"],
        "timing_base_domain_sha256": timing_domain["base_domain_sha256"],
        "timing_domain_sha256": timing["domain_sha256"],
        "timing_qualification_map_sha256":
            timing["timing_qualification_map_sha256"],
        "recovery_freeze_manifest_sha256": recovery_freeze,
        "timing_freeze_manifest_sha256": timing_freeze,
        "architecture_artifact_sha256": recovery_artifacts,
        "recovery_cells_per_arm": recovery_domain["expected_cells_per_arm"],
        "timing_rows": expected_timing_rows,
        "candidate_roster": candidates,
    }
    if not eligible:
        result = {
            "schema": SCHEMA + ".architecture-selection.v1",
            **bindings,
            "eligible_candidates": [],
            "eligible_overhead0_failures": {},
            "minimum_overhead0_failures": None,
            "recovery_equivalence_allowance": None,
            "recovery_equivalent_candidates": [],
            "selected_arm": None,
            "selected_codec": None,
            "selected_arm_descriptor_sha256": None,
            "selected_architecture_sha256": None,
            "ranking": [],
        }
        result["selection_sha256"] = sha256_json(result)
        return result
    minimum = min(failures[arm]["0"] for arm in eligible)
    margin_ppm = contract["selection"][
        "architecture_failure_equivalence_ppm"]
    cells = recovery_domain["expected_cells_per_arm"]
    allowance = max(1, (margin_ppm * cells + 1000000 - 1) // 1000000)
    equivalent = [
        arm for arm in eligible
        if failures[arm]["0"] <= minimum + allowance
    ]
    ranking_values = []
    for arm in equivalent:
        try:
            solve = timing_candidates[arm]["panels"][
                "decoder_solve|wirehair2_head"]["aggregate"]
        except (KeyError, TypeError):
            fail("development timing summary lacks aggregate solve evidence")
        if not isinstance(solve, dict):
            fail("development aggregate solve evidence is malformed")
        mean = solve.get("mean_log_ratio")
        if (solve.get("selectable") is not True or type(mean) is not float or
                not math.isfinite(mean)):
            fail("development aggregate solve evidence is not selectable")
        ranking_values.append((
            mean, failures[arm]["1"], failures[arm]["2"],
            failures[arm]["4"], arm,
        ))
    ranking_values.sort()
    selected = ranking_values[0][4]
    selected_artifact = arm_artifacts[selected]
    result = {
        "schema": SCHEMA + ".architecture-selection.v1",
        **bindings,
        "eligible_candidates": sorted(eligible),
        "eligible_overhead0_failures": {
            arm: failures[arm]["0"] for arm in sorted(eligible)
        },
        "minimum_overhead0_failures": minimum,
        "recovery_equivalence_allowance": allowance,
        "recovery_equivalent_candidates": sorted(equivalent),
        "ranking": [
            {
                "arm": value[4], "decoder_solve_mean_log_ratio": value[0],
                "failures_overhead0": failures[value[4]]["0"],
                "failures_overhead1": value[1],
                "failures_overhead2": value[2],
                "failures_overhead4": value[3],
            }
            for value in ranking_values
        ],
        "selected_arm": selected,
        "selected_codec": selected_artifact["codec"],
        "selected_arm_descriptor_sha256":
            selected_artifact["arm_descriptor_sha256"],
        "selected_architecture_sha256": selected_architecture_sha256(
            selected, selected_artifact),
    }
    result["selection_sha256"] = sha256_json(result)
    return result


def _histogram(values: Iterable[int]) -> Dict[str, int]:
    result: Dict[str, int] = defaultdict(int)
    for value in values:
        result[str(value)] += 1
    return dict(sorted(result.items(), key=lambda item: int(item[0])))


def _command_hashes(contract_path: Path) -> int:
    contract = load_contract(contract_path, check_domain_hashes=False)
    for phase in sorted(contract["recovery"]["domains"]):
        print("recovery:{} {}".format(
            phase, recovery_domain_sha256(contract, phase)))
    for phase in sorted(contract["timing"]["domains"]):
        print("timing-base:{} {}".format(
            phase, timing_base_domain_sha256(contract, phase)))
    return 0


def _command_describe(contract_path: Path) -> int:
    contract = load_contract(contract_path)
    result = {
        "schema": contract["schema"],
        "contract_id": contract["contract_id"],
        "contract_sha256": contract_sha256(contract),
        "recovery_domains": {
            phase: {
                "cells_per_arm": domain["expected_cells_per_arm"],
                "domain_sha256": domain["domain_sha256"],
            }
            for phase, domain in sorted(contract["recovery"]["domains"].items())
        },
        "timing_domains": {
            phase: {
                "cells": domain["expected_cells"],
                "base_domain_sha256": domain["base_domain_sha256"],
            }
            for phase, domain in sorted(contract["timing"]["domains"].items())
        },
        "development_wall_time_seconds_per_candidate":
            contract["timing"]["development_wall_time_seconds_per_candidate"],
    }
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


def _repair_map_arguments(values: Sequence[str]) -> Mapping[str, Path]:
    result: Dict[str, Path] = {}
    for value in values:
        arm, separator, path_text = value.partition("=")
        if not separator or not arm or not path_text or arm in result:
            fail("--repair-map must be a unique ARM=PATH pair")
        result[arm] = Path(path_text)
    return result


def _final_freeze_arguments(
        values: Sequence[str]) -> Mapping[Tuple[str, str], Path]:
    result: Dict[Tuple[str, str], Path] = {}
    for value in values:
        key_text, separator, path_text = value.partition("=")
        evidence_kind, kind_separator, phase = key_text.partition(":")
        key = (evidence_kind, phase)
        if (not separator or not kind_separator or not path_text or
                key in result):
            fail("--freeze must be a unique KIND:PHASE=PATH value")
        result[key] = Path(path_text)
    return result


def main(argv: Sequence[str] = ()) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("hashes", help="derive domain hashes while authoring a contract")
    subparsers.add_parser("describe", help="validate and summarize the frozen contract")
    ledger = subparsers.add_parser("validate-ledger", help="validate and summarize JSONL results")
    ledger.add_argument("--phase", required=True)
    ledger.add_argument("--freeze-manifest", required=True, type=Path)
    ledger.add_argument("--trace-manifest", required=True, type=Path)
    ledger.add_argument(
        "--repair-map", action="append", default=[], metavar="ARM=PATH",
        help="canonical frozen repair map; repeat once per mapped arm")
    ledger.add_argument("ledger", type=Path)
    timing = subparsers.add_parser(
        "validate-timing", help="validate and summarize paired timing JSONL")
    timing.add_argument("--phase", required=True)
    timing.add_argument("--freeze-manifest", required=True, type=Path)
    timing.add_argument("--trace-manifest", required=True, type=Path)
    timing.add_argument(
        "--timing-qualification-map", required=True, type=Path,
        help="canonical phase-scoped candidate-blind loss-retry map")
    timing.add_argument(
        "--timing-qualification-audit", required=True, type=Path,
        help="ordered canonical audit of every attempted loss retry")
    timing.add_argument(
        "--timing-qualification-map-sha256", required=True,
        help="pre-result content identity of the qualification map")
    timing.add_argument(
        "--repair-map", action="append", default=[], metavar="ARM=PATH",
        help="canonical frozen repair map; repeat once per mapped arm")
    timing.add_argument("receipt", type=Path)
    continuity = subparsers.add_parser(
        "validate-final-continuity",
        help="prove that every final phase uses one architecture and map set")
    continuity.add_argument(
        "--freeze", action="append", default=[],
        metavar="KIND:PHASE=PATH",
        help="repeat for all four final recovery freezes and final timing")
    continuity.add_argument(
        "--selection-receipt", required=True, type=Path,
        help="canonical development architecture-selection receipt")
    continuity.add_argument(
        "--timing-qualification-map", required=True, type=Path)
    continuity.add_argument(
        "--timing-qualification-audit", required=True, type=Path)
    continuity.add_argument(
        "--timing-qualification-map-sha256", required=True)
    continuity.add_argument(
        "--timing-trace-manifest", required=True, type=Path,
        help="final ordered timing trace manifest frozen before results")
    args = parser.parse_args(list(argv) if argv else None)
    try:
        if args.command == "hashes":
            return _command_hashes(args.contract)
        if args.command == "describe":
            return _command_describe(args.contract)
        contract = load_contract(args.contract)
        if args.command == "validate-final-continuity":
            timing_qualification = load_timing_qualification_map(
                contract, "final", args.timing_qualification_map,
                args.timing_qualification_audit,
                args.timing_qualification_map_sha256)
            summary = validate_final_freeze_continuity(
                contract, _final_freeze_arguments(args.freeze),
                _load_canonical_json_file(
                    args.selection_receipt, "architecture selection receipt"),
                timing_qualification, args.timing_trace_manifest)
            print(json.dumps(summary, sort_keys=True, indent=2))
            return 0
        repair_maps = _repair_map_arguments(args.repair_map)
        if args.command == "validate-ledger":
            summary = validate_ledger(
                contract, args.phase, args.ledger, args.freeze_manifest,
                args.trace_manifest, repair_maps)
        else:
            summary = validate_timing_receipt(
                contract, args.phase, args.receipt, args.freeze_manifest,
                args.trace_manifest, args.timing_qualification_map,
                args.timing_qualification_audit,
                args.timing_qualification_map_sha256, repair_maps)
        print(json.dumps(summary, sort_keys=True, indent=2))
        return 0
    except ContractError as exc:
        print("contract error: {}".format(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
