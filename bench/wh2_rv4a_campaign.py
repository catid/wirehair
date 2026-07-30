#!/usr/bin/env python3
"""Authenticated result-free campaign for WH2 repair-v1 (bead rv4a.2).

The training and joint sealed populations in this module are immutable and
are constructed without reading a codec result.  ``plan`` is therefore safe
before a benchmark binary exists.  Execution and verification entry points
are installed only with the complete v3 parser and otherwise fail closed.

The joint sealed phase deliberately has no lane or survivor command-line
switch.  It consumes the sole eligible winner recorded by a completed
training campaign and executes both prebound sealed lanes into one completion.
"""

import argparse
from contextlib import contextmanager
import functools
import gzip
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path
import queue
import signal
import stat
import subprocess
import sys
import tempfile
import time
import types
import zlib


REPOSITORY = Path(__file__).resolve().parents[1]
TOOLS = REPOSITORY / "tools"
_BOOTSTRAP_SOURCE_PINS = globals().get("_BOOTSTRAP_SOURCE_PINS")
_BOOTSTRAP_RUNNER_SOURCE_SHA256 = globals().get(
    "_BOOTSTRAP_RUNNER_SOURCE_SHA256")
# Never ask the import machinery to choose between source and a timestamp-
# valid bytecode cache for the two local evidence consumers.  Compile stable,
# no-follow source snapshots directly.  This remains safe even if an attacker
# pre-populates both __pycache__ and a deterministic alternate cache tree.
sys.dont_write_bytecode = True
sys.path[:] = [entry for entry in sys.path if entry != str(TOOLS)]
sys.path.insert(0, str(TOOLS))


def _bootstrap_source_bytes(path, *, byte_limit=2 * 1024 * 1024):
    """Snapshot one local Python source without following or reopening it."""
    path = Path(path).resolve()
    flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NONBLOCK", 0)
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = None
    try:
        descriptor = os.open(str(path), flags)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode) or
            before.st_size > byte_limit
        ):
            raise RuntimeError(f"invalid local Python source: {path}")
        chunks = []
        size = 0
        while True:
            chunk = os.read(
                descriptor, min(64 * 1024, byte_limit - size + 1))
            if not chunk:
                break
            size += len(chunk)
            if size > byte_limit:
                raise RuntimeError(f"oversized local Python source: {path}")
            chunks.append(chunk)
        after = os.fstat(descriptor)
        current = os.stat(path, follow_symlinks=False)
    finally:
        if descriptor is not None:
            os.close(descriptor)
    stable = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    if (
        any(getattr(before, name) != getattr(after, name)
            for name in stable) or
        any(getattr(after, name) != getattr(current, name)
            for name in stable) or
        size != after.st_size
    ):
        raise RuntimeError(f"unstable local Python source: {path}")
    return b"".join(chunks)


def _bootstrap_source_module(module_name, path):
    """Execute exact source bytes under the public dependency module name."""
    path = Path(path).resolve()
    payload = _bootstrap_source_bytes(path)
    payload_sha256 = hashlib.sha256(payload).hexdigest()
    if _BOOTSTRAP_SOURCE_PINS is not None:
        if (
            not isinstance(_BOOTSTRAP_SOURCE_PINS, dict) or
            set(_BOOTSTRAP_SOURCE_PINS) != {
                "peel_codec", "repair_timing_codec"
            } or
            _BOOTSTRAP_SOURCE_PINS.get(module_name) != payload_sha256
        ):
            raise RuntimeError(
                f"local Python source does not match bootstrap pin: {path}")
    module = types.ModuleType(module_name)
    module.__file__ = str(path)
    module.__package__ = ""
    module.__loader__ = None
    module.__spec__ = None
    missing = object()
    previous = sys.modules.get(module_name, missing)
    sys.modules[module_name] = module
    try:
        code = compile(payload, str(path), "exec", dont_inherit=True)
        exec(code, module.__dict__)
        module.__rv4a_source_sha256__ = payload_sha256
        if _bootstrap_source_bytes(path) != payload:
            raise RuntimeError(
                f"local Python source changed during import: {path}")
    except BaseException:
        if previous is missing:
            sys.modules.pop(module_name, None)
        else:
            sys.modules[module_name] = previous
        raise
    return module


peel_codec = _bootstrap_source_module(
    "peel_codec", TOOLS / "peel_codec.py")

# Keep result-free planning usable while the additive parser and native pieces
# are being built in parallel.  Execution and verification fail closed unless
# the exact v3 parser source is present.
try:
    repair = _bootstrap_source_module(
        "repair_timing_codec", TOOLS / "repair_timing_codec.py")
except FileNotFoundError:  # pragma: no cover - pre-integration only
    repair = None


CAMPAIGN_SCHEMA = "wirehair.wh2.rv4a.repair-campaign.v1"
JOB_SCHEMA = "wirehair.wh2.rv4a.repair-job.v1"
RESULT_SCHEMA = "wirehair.wh2.rv4a.repair-result.v1"
SUMMARY_SCHEMA = "wirehair.wh2.rv4a.repair-summary.v1"
COMPLETION_SCHEMA = "wirehair.wh2.rv4a.repair-completion.v1"
TRAINING_DECISION_SCHEMA = "wirehair.wh2.rv4a.training-decision.v1"
REQUIRED_REPAIRTIMING_PROTOCOL = (
    "wirehair-v2-bench:repairtiming:repair-v1:v3"
)
REQUIRED_REPAIRTIMING_SCHEMA = "wirehair.wh2.repairtiming.v3"

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
CORRECTNESS_BLOCK_BYTES = (2, 6, 32, 64, 256, 1280, 4096)
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
REPAIR_ATTEMPTS = 8
SELECTOR_COUNT_FIELDS = (
    "attempts_executed",
    "calls_executed",
    "real_configuration_calls",
    "structural_probe_calls",
)
SELECTOR_COUNT_LIMITS = {
    "attempts_executed": (1, REPAIR_ATTEMPTS),
    "calls_executed": (1, REPAIR_ATTEMPTS + 2),
    "real_configuration_calls": (1, 2),
    "structural_probe_calls": (0, REPAIR_ATTEMPTS),
}
TRAINING_REPLICATES = 768
SEALED_REPLICATES = 256

CONSTRUCTION_DERIVED = (
    "base_xor_d192ed03_times_rep_plus_1_mod2^32_v1"
)
CONSTRUCTION_FIXED = "fixed_base_u32_v1"
LOSS_DERIVED = (
    "base_xor_9e3779b97f4a7c15_times_rep_plus_1_v1"
)
PRECODE_ATTEMPT_INCREMENT = 0x9E3779B97F4A7C15
PACKET_ATTEMPT_INCREMENT = 0x9E3779B9

REPAIR_POLICY_SHA256 = (
    "5e67a150d1f909d6ed80468185fa2dd0"
    "e82eb2fc3486c0fa662e213cf3100b42"
)
ARMS = (
    {
        "name": "pure8_s0m1_d3",
        "native_name": "pure8_s0m1_d3_repair_v1",
        "provisional_id": "19cccf775ce0bf09",
        "contract_sha256":
            "19cccf775ce0bf098c9a425cb349714c"
            "4c4a880e7cf136c3bc365e13c05089a5",
        "gf256_rows": 8,
        "gf16_rows": 0,
        "dense_rows": 3,
        "staircase": "SmallBandStaircaseCount(K)-1",
        "period": 244,
        "x_mode": "frozen",
        "packet_mix_count": 3,
        "historical_za5v_raw_weak_constructions": 1208,
    },
    {
        "name": "pure9_s0m1_d3",
        "native_name": "pure9_s0m1_d3_repair_v1",
        "provisional_id": "a530f9105beaa450",
        "contract_sha256":
            "a530f9105beaa450dee70ad9b2a5cc54"
            "c944d3cd47f0aa6534630b8971608541",
        "gf256_rows": 9,
        "gf16_rows": 0,
        "dense_rows": 3,
        "staircase": "SmallBandStaircaseCount(K)-1",
        "period": 244,
        "x_mode": "frozen",
        "packet_mix_count": 3,
        "historical_za5v_raw_weak_constructions": 1034,
    },
)
ARM_BY_NAME = {arm["name"]: arm for arm in ARMS}
ARM_BY_ID = {arm["provisional_id"]: arm for arm in ARMS}

# Prior result-free populations whose roots/loss seeds may never be reused.
ZA5V_SELECTION_CONSTRUCTION_BASE = 0x13579BDF
ZA5V_SELECTION_LOSS_BASE = 0x0123456789ABCDEF
ZA5V_HOLDOUT_CONSTRUCTION_BASE = 0x2468ACE0
ZA5V_HOLDOUT_LOSS_BASE = 0xFEDCBA9876543210
ZA5V_SCREEN_CONSTRUCTION_BASE = 0x6789CAFE
ZA5V_SCREEN_LOSS_BASE = 0x243F6A8885A308D3

SEED_KDF_SCHEMA = "wirehair.wh2.rv4a.domain-kdf.v1"
SEED_KDF_PREFIX = b"wirehair.wh2.rv4a.domain-kdf.v1\x00"
SEED_DOMAINS = {
    "training_construction":
        "wirehair.wh2.rv4a.training.construction-base.v1",
    "sealed_random_construction":
        "wirehair.wh2.rv4a.sealed-random.construction-base.v1",
    "sealed_production_construction":
        "wirehair.wh2.rv4a.sealed-production.construction-root.v1",
    "training_loss":
        "wirehair.wh2.rv4a.training.loss-base.v1",
    "sealed_random_loss":
        "wirehair.wh2.rv4a.sealed-random.loss-base.v1",
    "sealed_production_loss":
        "wirehair.wh2.rv4a.sealed-production.loss-base.v1",
    "training_order": "wirehair.wh2.rv4a.training.job-order.v1",
    "sealed_order": "wirehair.wh2.rv4a.sealed.job-order.v1",
}

# Result envelopes contain one bounded normalized v3 receipt.  The exact
# parser-derived per-cell ceiling is checked at runtime; these fallback
# constants keep malformed evidence bounded even if a wrong parser is loaded.
# Once the parser is integrated, validate_frozen_contract() requires its
# ceiling to fit within this explicit campaign cap.
RESULT_UNCOMPRESSED_FIXED_BYTES = 2 * 1024 * 1024
RESULT_UNCOMPRESSED_BYTES_PER_CELL = 384 * 1024
SUMMARY_COMPRESSED_BYTE_LIMIT = 64 * 1024 * 1024
SUMMARY_UNCOMPRESSED_BYTE_LIMIT = 512 * 1024 * 1024
BENCHMARK_BINARY_BYTE_LIMIT = 512 * 1024 * 1024
NATIVE_SOURCE_BYTE_LIMIT = 8 * 1024 * 1024
PYTHON_SOURCE_BYTE_LIMIT = 2 * 1024 * 1024
BUILD_SOURCE_BYTE_LIMIT = 4 * 1024 * 1024
BUILD_CONFIGURATION_BYTE_LIMIT = 16 * 1024 * 1024
RUNTIME_PIN_SCHEMA = "wirehair.wh2.rv4a.verifier-pins.v1"
RUNTIME_PIN_BYTE_LIMIT = 2 * 1024 * 1024
REPLAY_QUEUE_MAX_ITEMS = 2 * MAX_WORKERS
REPLAY_QUEUE_POLL_SECONDS = 0.05
REPLAY_TERMINATE_GRACE_SECONDS = 0.5
REPLAY_KILL_GRACE_SECONDS = 1.0
REPLAY_SHUTDOWN_GRACE_SECONDS = 5.0
WORKER_POLL_SECONDS = 0.05
WORKER_TERMINATE_GRACE_SECONDS = 5.0
WORKER_KILL_GRACE_SECONDS = 1.0

_PARENT_SIGNAL_STATE = None


class CampaignError(RuntimeError):
    """The frozen campaign contract or authenticated evidence was invalid."""


def canonical_json_bytes(value):
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


def _is_sha256(value):
    return (
        type(value) is str and len(value) == 64 and
        all(character in "0123456789abcdef" for character in value)
    )


def construction_seed(base, replicate):
    if (
        type(base) is not int or not 0 <= base <= 0xFFFFFFFF or
        type(replicate) is not int or replicate < 0
    ):
        raise CampaignError("invalid derived construction seed input")
    return (
        base ^ (((replicate + 1) * 0xD192ED03) & 0xFFFFFFFF)
    ) & 0xFFFFFFFF


def loss_seed(base, replicate):
    if (
        type(base) is not int or not 0 <= base <= 0xFFFFFFFFFFFFFFFF or
        type(replicate) is not int or replicate < 0
    ):
        raise CampaignError("invalid derived loss seed input")
    return (
        base ^
        (((replicate + 1) * 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF)
    ) & 0xFFFFFFFFFFFFFFFF


def expanded_precode_seed(root, attempt):
    return (
        root + attempt * PRECODE_ATTEMPT_INCREMENT
    ) & 0xFFFFFFFFFFFFFFFF


def expanded_packet_seed(root, attempt):
    return (root + attempt * PACKET_ATTEMPT_INCREMENT) & 0xFFFFFFFF


def _kdf_payload(domain, coordinates, rejection_counter):
    if (
        type(domain) is not str or not domain.isascii() or not domain or
        type(coordinates) not in (tuple, list) or
        any(type(value) is not int or not 0 <= value <= 0xFFFFFFFFFFFFFFFF
            for value in coordinates) or
        type(rejection_counter) is not int or
        not 0 <= rejection_counter <= 0xFFFFFFFF
    ):
        raise CampaignError("invalid seed-domain KDF input")
    encoded = domain.encode("ascii")
    if len(encoded) > 0xFFFF or len(coordinates) > 0xFFFF:
        raise CampaignError("seed-domain KDF framing overflow")
    return b"".join((
        SEED_KDF_PREFIX,
        len(encoded).to_bytes(2, "big"),
        encoded,
        len(coordinates).to_bytes(2, "big"),
        b"".join(value.to_bytes(8, "big") for value in coordinates),
        rejection_counter.to_bytes(4, "big"),
    ))


def domain_value(domain, width_bits, *, coordinates=(), counter=0):
    if width_bits not in (32, 64):
        raise CampaignError("seed-domain KDF width must be 32 or 64")
    digest = hashlib.sha256(
        _kdf_payload(domain, coordinates, counter)).digest()
    return int.from_bytes(digest[:width_bits // 8], "big")


def _expanded_sets(roots):
    precode = set()
    packet = set()
    for root in roots:
        for attempt in range(REPAIR_ATTEMPTS):
            precode.add(expanded_precode_seed(root, attempt))
            packet.add(expanded_packet_seed(root, attempt))
    return precode, packet


def _root_population(base, count, derivation):
    if derivation == CONSTRUCTION_DERIVED:
        return tuple(construction_seed(base, index) for index in range(count))
    if derivation == CONSTRUCTION_FIXED:
        return tuple(base for unused_index in range(count))
    raise CampaignError("unknown construction derivation")


def _loss_population(base, count):
    return tuple(loss_seed(base, index) for index in range(count))


def _digest_ints(values):
    return canonical_sha256(list(values))


def _seed_domain_contract():
    """Derive and exhaustively prove every pre-outcome seed population.

    Domain candidates are deterministic.  A whole lane (or one production K)
    is rejected and its counter incremented if either expanded attempt seed
    space collides with a prior or already accepted lane.  Roots intentionally
    reused across arms, K values, and schedules are represented once here.
    """
    prior_roots = {
        "za5v_selection": _root_population(
            ZA5V_SELECTION_CONSTRUCTION_BASE, 768, CONSTRUCTION_DERIVED),
        "za5v_forbidden_holdout": _root_population(
            ZA5V_HOLDOUT_CONSTRUCTION_BASE, 256, CONSTRUCTION_DERIVED),
        "za5v_frozen_and_neighbor_screens": _root_population(
            ZA5V_SCREEN_CONSTRUCTION_BASE, 3, CONSTRUCTION_DERIVED),
    }
    prior_losses = {
        "za5v_selection": _loss_population(
            ZA5V_SELECTION_LOSS_BASE, 768),
        "za5v_forbidden_holdout": _loss_population(
            ZA5V_HOLDOUT_LOSS_BASE, 256),
        "za5v_frozen_and_neighbor_screens": _loss_population(
            ZA5V_SCREEN_LOSS_BASE, 3),
    }

    # Historical campaigns executed attempt zero only.  Their actual seed
    # evidence is therefore the zero-extended root in both seed spaces, not a
    # fictitious set of attempts 1..7.  New lanes expand every root through
    # all eight repair-v1 attempts and are rejected against these actual sets.
    occupied_precode = set()
    occupied_packet = set()
    conservative_prior_precode = set()
    conservative_prior_packet = set()
    prior_root_receipts = {}
    for name, roots in prior_roots.items():
        if len(set(roots)) != len(roots):
            raise CampaignError(f"prior root population is not unique: {name}")
        actual_precode = {expanded_precode_seed(root, 0) for root in roots}
        actual_packet = {expanded_packet_seed(root, 0) for root in roots}
        precode, packet = _expanded_sets(roots)
        if (
            len(precode) != len(roots) * REPAIR_ATTEMPTS or
            len(packet) != len(roots) * REPAIR_ATTEMPTS
        ):
            raise CampaignError(
                f"prior expanded attempt population collides: {name}")
        if (
            occupied_precode & actual_precode or
            occupied_packet & actual_packet
        ):
            raise CampaignError("prior actual attempt-zero populations overlap")
        occupied_precode.update(actual_precode)
        occupied_packet.update(actual_packet)
        conservative_prior_precode.update(precode)
        conservative_prior_packet.update(packet)
        prior_root_receipts[name] = {
            "count": len(roots),
            "root_set_sha256": _digest_ints(sorted(roots)),
            "actual_attempt0_precode64_set_sha256":
                _digest_ints(sorted(actual_precode)),
            "actual_attempt0_packet32_set_sha256":
                _digest_ints(sorted(actual_packet)),
            "conservative_hypothetical_attempts0_through7_precode64_sha256":
                _digest_ints(sorted(precode)),
            "conservative_hypothetical_attempts0_through7_packet32_sha256":
                _digest_ints(sorted(packet)),
        }

    accepted_roots = {}

    def accept_derived_lane(name, count):
        domain = SEED_DOMAINS[name]
        for counter in range(0x100000000):
            base = domain_value(domain, 32, counter=counter)
            roots = _root_population(base, count, CONSTRUCTION_DERIVED)
            if len(set(roots)) != count:
                continue
            precode, packet = _expanded_sets(roots)
            if (
                len(precode) != count * REPAIR_ATTEMPTS or
                len(packet) != count * REPAIR_ATTEMPTS or
                occupied_precode & precode or
                occupied_packet & packet
            ):
                continue
            occupied_precode.update(precode)
            occupied_packet.update(packet)
            accepted_roots[name] = {
                "domain": domain,
                "coordinates": [],
                "rejection_count": counter,
                "base": base,
                "derivation": CONSTRUCTION_DERIVED,
                "count": count,
                "root_set_sha256": _digest_ints(sorted(roots)),
                "ordered_root_list_sha256": _digest_ints(roots),
                "expanded_precode64_set_sha256":
                    _digest_ints(sorted(precode)),
                "expanded_packet32_set_sha256":
                    _digest_ints(sorted(packet)),
                "_roots": roots,
            }
            return
        raise CampaignError(f"could not derive collision-free lane {name}")

    accept_derived_lane("training_construction", TRAINING_REPLICATES)
    accept_derived_lane(
        "sealed_random_construction", SEALED_REPLICATES)

    production = {}
    production_domain = SEED_DOMAINS["sealed_production_construction"]
    for block_count in K_VALUES:
        for counter in range(0x100000000):
            root = domain_value(
                production_domain,
                32,
                coordinates=(block_count,),
                counter=counter,
            )
            precode, packet = _expanded_sets((root,))
            if (
                len(precode) != REPAIR_ATTEMPTS or
                len(packet) != REPAIR_ATTEMPTS or
                occupied_precode & precode or
                occupied_packet & packet
            ):
                continue
            occupied_precode.update(precode)
            occupied_packet.update(packet)
            production[str(block_count)] = {
                "domain": production_domain,
                "coordinates": [block_count],
                "rejection_count": counter,
                "root": root,
                "derivation": CONSTRUCTION_FIXED,
                "count": 1,
                "expanded_precode64_count": len(precode),
                "expanded_packet32_count": len(packet),
                "expanded_precode64_set_sha256":
                    _digest_ints(sorted(precode)),
                "expanded_packet32_set_sha256":
                    _digest_ints(sorted(packet)),
            }
            break
        else:
            raise CampaignError(
                f"could not derive production root for K={block_count}")
    production_roots = tuple(
        production[str(block_count)]["root"] for block_count in K_VALUES)
    if len(set(production_roots)) != len(production_roots):
        raise CampaignError("production roots are not unique")

    occupied_losses = set()
    prior_loss_receipts = {}
    for name, values in prior_losses.items():
        if len(set(values)) != len(values) or occupied_losses & set(values):
            raise CampaignError("prior loss populations overlap")
        occupied_losses.update(values)
        prior_loss_receipts[name] = {
            "count": len(values),
            "loss_set_sha256": _digest_ints(sorted(values)),
        }

    accepted_losses = {}
    for name, count in (
        ("training_loss", TRAINING_REPLICATES),
        ("sealed_random_loss", SEALED_REPLICATES),
        ("sealed_production_loss", SEALED_REPLICATES),
    ):
        domain = SEED_DOMAINS[name]
        for counter in range(0x100000000):
            base = domain_value(domain, 64, counter=counter)
            values = _loss_population(base, count)
            if len(set(values)) != count or occupied_losses & set(values):
                continue
            occupied_losses.update(values)
            accepted_losses[name] = {
                "domain": domain,
                "coordinates": [],
                "rejection_count": counter,
                "base": base,
                "derivation": LOSS_DERIVED,
                "count": count,
                "loss_set_sha256": _digest_ints(sorted(values)),
                "ordered_loss_list_sha256": _digest_ints(values),
                "_losses": values,
            }
            break
        else:
            raise CampaignError(
                f"could not derive collision-free loss lane {name}")

    orders = {
        phase: {
            "domain": SEED_DOMAINS[f"{phase}_order"],
            "coordinates": [],
            "rejection_count": 0,
            "value": domain_value(
                SEED_DOMAINS[f"{phase}_order"], 64),
        }
        for phase in ("training", "sealed")
    }
    if orders["training"]["value"] == orders["sealed"]["value"]:
        raise CampaignError("job-order domains collided")

    def public_receipt(item):
        return {
            key: value for key, value in item.items()
            if not key.startswith("_")
        }

    return {
        "schema": "wirehair.wh2.rv4a.seed-disjointness.v1",
        "kdf": {
            "schema": SEED_KDF_SCHEMA,
            "hash": "sha256",
            "framing":
                "prefix+u16be-domain-bytes+u16be-coordinate-count+"
                "u64be-coordinates+u32be-rejection-counter",
            "output": "most-significant-32-or-64-bits",
            "rejection":
                "increment-u32be-counter-until-entire-lane-has-unique-"
                "roots-losses-and-separate-expanded-precode64-packet32-"
                "sets-disjoint-from-prior-actual-attempt0-and-earlier-"
                "new-lanes",
        },
        "construction": {
            "prior": prior_root_receipts,
            "training": public_receipt(
                accepted_roots["training_construction"]),
            "sealed_random": public_receipt(
                accepted_roots["sealed_random_construction"]),
            "sealed_production": {
                "domain": production_domain,
                "derivation": CONSTRUCTION_FIXED,
                "roots_by_K": production,
                "root_count": len(production_roots),
                "root_set_sha256":
                    _digest_ints(sorted(production_roots)),
                "ordered_K_root_list_sha256":
                    _digest_ints(production_roots),
            },
            "new_expanded_precode64_vs_prior_actual_and_new_intersections": 0,
            "new_expanded_packet32_vs_prior_actual_and_new_intersections": 0,
            "prior_evidence_scope": "actual-attempt0-only",
            "conservative_hypothetical_prior_expansion": {
                "not_historical_evidence": True,
                "new_precode64_intersection_count": len(
                    conservative_prior_precode &
                    (occupied_precode - {
                        expanded_precode_seed(root, 0)
                        for roots in prior_roots.values() for root in roots
                    })
                ),
                "new_packet32_intersection_count": len(
                    conservative_prior_packet &
                    (occupied_packet - {
                        expanded_packet_seed(root, 0)
                        for roots in prior_roots.values() for root in roots
                    })
                ),
            },
        },
        "loss": {
            "prior": prior_loss_receipts,
            "training": public_receipt(
                accepted_losses["training_loss"]),
            "sealed_random": public_receipt(
                accepted_losses["sealed_random_loss"]),
            "sealed_production": public_receipt(
                accepted_losses["sealed_production_loss"]),
            "all_cross_population_intersections": 0,
        },
        "order": orders,
        "_private": {
            "training_roots":
                accepted_roots["training_construction"]["_roots"],
            "sealed_random_roots":
                accepted_roots["sealed_random_construction"]["_roots"],
            "production_roots": production_roots,
            "training_losses":
                accepted_losses["training_loss"]["_losses"],
            "sealed_random_losses":
                accepted_losses["sealed_random_loss"]["_losses"],
            "sealed_production_losses":
                accepted_losses["sealed_production_loss"]["_losses"],
        },
    }


_SEEDS = _seed_domain_contract()


def seed_disjointness_proof():
    return {
        key: value for key, value in _SEEDS.items()
        if key != "_private"
    }


def selection_policy():
    """Return the prebound policy; this function never reads outcomes."""
    # "Available" is relative to the frozen repairtiming v3 stats tuple, not
    # an instruction to reinterpret that schema when PrecodeSolveStats grows.
    work_statistics = {
        "fields":
            "every-available-PrecodeSolveStats-work-counter-and-phase-ns",
        "availability":
            "aggregate-available-values-only-and-report-available-and-"
            "unavailable-counts-never-fold-unavailable-into-zero",
        "statistics": ["mean", "p50", "p95", "p99", "max"],
        "percentile":
            "nearest-rank-on-nondecreasing-integer-values-rank-ceil(p*n)-1",
    }
    timing_confirmation = {
        "job_metric":
            "arithmetic-mean-of-eligible-paired-replicate-log-costs",
        "aggregate":
            "equal-weight-arithmetic-mean-over-all-K-schedule-jobs",
        "uncertainty":
            "max(two-sided-Student-t-95-percent-upper-bound-over-job-"
            "means,mean-of-authenticated-per-job-CI-highs)",
        "pass":
            "aggregate-upper-strictly-less-than-negative-mean-effective-"
            "floor-including-required-margin-and-both-corresponding-AA-"
            "noise-floors",
        "completeness":
            "every-job-cross-contrast-and-both-corresponding-AA-arms-"
            "finite-with-at-least-half-the-frozen-replicates",
    }
    return {
        "schema": "wirehair.wh2.rv4a.repair-selection-policy.v1",
        "scope": {
            "training":
                "both-frozen-arms-all-K2-through100-all-six-schedules-"
                "768-paired-roots-block-bytes2",
            "sealed":
                "one-eligible-training-winner-only-two-joint-lanes-all-"
                "K2-through100-all-six-schedules-256-loss-replicates-"
                "block-bytes2-one-completion-no-fallback",
            "wider_block_bytes":
                "correctness-prelaunch-only-never-roster-or-ranking",
        },
        "raw_evidence": {
            "attempt_zero":
                "fresh-attempt0-result-and-recovery-metrics-reported-"
                "separately-with-no-seed-fixes-including-deduplicated-"
                "raw-and-repaired-weak-construction-counts",
            "historical_references": {
                arm["name"]:
                    arm["historical_za5v_raw_weak_constructions"]
                for arm in ARMS
            },
            "rule":
                "historical-counts-are-reference-only-fresh-v3-evidence-"
                "is-authoritative-and-is-never-erased-by-repaired-metrics",
        },
        "selector": {
            "policy_sha256": REPAIR_POLICY_SHA256,
            "attempts": list(range(REPAIR_ATTEMPTS)),
            "deduplication_key": ["arm_provisional_id", "K", "root"],
            "structural_duplicate":
                "exact-selector-projection-only-excludes-elapsed-payload-"
                "and-attempt0-real-stats",
            "schedule_rule":
                "same-key-selector-projection-must-be-identical-across-"
                "all-six-schedules",
            "weighting":
                "exactly-once-per-deduplication-key-never-per-physical-"
                "recovery-cell",
            "real_work_witness":
                "lowest-replicate-carrying-the-key-in-the-first-frozen-"
                "schedule",
            "execution_counts":
                "attempts-calls-real-configurations-and-structural-probes-"
                "each-reported-as-histogram-mean-p50-p95-p99-max-over-"
                "deduplicated-selectors",
            "weak_constructions":
                "fresh-raw-attempt0-and-repaired-final-weak-counts-once-per-"
                "deduplication-key",
            "work": work_statistics,
            "selected_attempt_histogram":
                "counts-attempt0-through7-plus-no-attempt-with-mean-p50-"
                "p95-p99-max-over-selected-attempt-values",
        },
        "recovery": {
            "matched_arms": ["candidate_repaired", "dispatch", "wh1"],
            "preserved_reference": "candidate_raw_attempt0",
            "final_unrecovered_h64":
                "sum-count-not-success-after-K-plus64-received-symbols",
            "tail_auc_h0_h64":
                "sum-over-integer-h0-through64-of-unrecovered-count",
            "weighting": "equal-weight-per-physical-recovery-cell",
            "nonrepairable_regression":
                "candidate-nonweak-failure-when-either-matched-dispatch-"
                "or-wh1-recovers-shared-failures-are-not-regressions",
        },
        "eligibility": {
            "zero": [
                "selector_cap_exhaustion",
                "selector_fatal_attempt_zero_mismatch",
                "selector_uncommitted",
                "candidate_final_weak",
                "candidate_nonrepairable_regression",
                "candidate_runtime_error",
                "selector_or-codec_OOM",
                "forced_control_failure",
            ],
            "recovery":
                "candidate-repaired-final-unrecovered-h64-and-tail-AUC-"
                "each-less-than-or-equal-to-both-matched-dispatch-and-wh1",
            "selected_direct_vs_dispatch": timing_confirmation,
            "full_encoder_vs_wh1": timing_confirmation,
            "decoder_feed_vs_wh1": timing_confirmation,
            "full_decoder_vs_wh1": timing_confirmation,
            "selector_vs_forced_encoder_reporting":
                "same-authenticated-log-cost-and-AA-noise-reporting-as-"
                "timing-confirmations-but-never-an-eligibility-gate",
            "direct_prefix":
                "pair-local-max-repaired-dispatch-first-success-when-both-"
                "recover-otherwise-K-plus64-shared-by-both-direct-solves",
            "no_eligible":
                "record-no-survivor-retain-dispatch-and-forbid-sealed",
        },
        "ranking": {
            "among": "eligible-training-arms-only",
            "lexicographic": [
                "decoder_feed_candidate_wh1_log_cost",
                "full_decoder_candidate_wh1_log_cost",
                "selected_direct_candidate_dispatch_log_cost",
                "full_encoder_candidate_wh1_log_cost",
                "candidate_name",
            ],
            "metric":
                "aggregate-equal-job-weight-log-cost-mean-smaller-is-faster",
        },
        "sealed": {
            "winner":
                "exactly-the-sole-training-decision-winner-no-CLI-choice",
            "lanes": [
                "256-disjoint-random-roots-with-256-paired-loss-seeds",
                "one-prebound-production-root-per-K-repeated-over-256-"
                "disjoint-loss-seeds",
            ],
            "confirmation":
                "repeat-all-training-eligibility-and-throughput-gates-on-"
                "the-union-and-independently-on-each-prebound-lane-in-one-"
                "joint-completion-with-no-adaptive-lane-choice",
            "failure": "retain-dispatch-no-fallback-arm-or-lane-inspection",
            "promotion":
                "rv4a2-records-evidence-only-public-contract-promotion-is-"
                "outside-this-task",
        },
        "thermal_edac": {
            "coverage":
                "every-worker-receipt-bound-to-assigned-singleton-CPU-and-"
                "manifest-thermal-device-inode-with-valid-samples-covering-"
                "the-native-work-interval",
            "errors":
                "corrected-and-uncorrected-EDAC-deltas-must-both-be-zero",
        },
    }


def _phase_parameters(phase, lane, block_count):
    construction = _SEEDS["construction"]
    losses = _SEEDS["loss"]
    if phase == "training" and lane == "training":
        return {
            "replicates": TRAINING_REPLICATES,
            "construction_seed_base": construction["training"]["base"],
            "construction_seed_derivation": CONSTRUCTION_DERIVED,
            "loss_seed_base": losses["training"]["base"],
            "loss_seed_derivation": LOSS_DERIVED,
        }
    if phase == "sealed" and lane == "random":
        return {
            "replicates": SEALED_REPLICATES,
            "construction_seed_base":
                construction["sealed_random"]["base"],
            "construction_seed_derivation": CONSTRUCTION_DERIVED,
            "loss_seed_base": losses["sealed_random"]["base"],
            "loss_seed_derivation": LOSS_DERIVED,
        }
    if phase == "sealed" and lane == "production":
        return {
            "replicates": SEALED_REPLICATES,
            "construction_seed_base":
                construction["sealed_production"]["roots_by_K"][
                    str(block_count)]["root"],
            "construction_seed_derivation": CONSTRUCTION_FIXED,
            "loss_seed_base": losses["sealed_production"]["base"],
            "loss_seed_derivation": LOSS_DERIVED,
        }
    raise CampaignError(f"invalid phase/lane {phase!r}/{lane!r}")


def _job_without_identity(phase, lane, arm_name, block_count, schedule):
    if (
        arm_name not in ARM_BY_NAME or block_count not in K_VALUES or
        schedule not in SCHEDULES
    ):
        raise CampaignError("job coordinate is outside frozen roster")
    arm = ARM_BY_NAME[arm_name]
    seeds = _phase_parameters(phase, lane, block_count)
    return {
        "schema": JOB_SCHEMA,
        "phase": phase,
        "lane": lane,
        "arm": arm_name,
        "arm_native_name": arm["native_name"],
        "arm_provisional_id": arm["provisional_id"],
        "arm_contract_sha256": arm["contract_sha256"],
        "repair_policy_sha256": REPAIR_POLICY_SHA256,
        "K": block_count,
        "schedule": schedule,
        "block_bytes": BLOCK_BYTES,
        "loss": LOSS,
        **seeds,
        "warmup_replicates": WARMUP_REPLICATES,
        "inner_reps": INNER_REPS,
        "max_overhead": MAX_OVERHEAD,
        "cache_state": CACHE_STATE,
        "systematic_cache": SYSTEMATIC_CACHE,
        "evict_bytes": EVICT_BYTES,
        "required_margin": REQUIRED_MARGIN,
        "thermal_sampling_interval_ms": THERMAL_INTERVAL_MS,
    }


def build_pre_cpu_jobs(phase, survivor=None):
    if phase == "training":
        if survivor is not None:
            raise CampaignError("training forbids a survivor")
        arm_names = tuple(arm["name"] for arm in ARMS)
        lanes = ("training",)
    elif phase == "sealed":
        if survivor not in ARM_BY_NAME:
            raise CampaignError("sealed requires one frozen survivor")
        arm_names = (survivor,)
        lanes = ("random", "production")
    else:
        raise CampaignError(f"unknown phase {phase!r}")

    jobs = []
    for arm_name in arm_names:
        for lane in lanes:
            for block_count in K_VALUES:
                for schedule in SCHEDULES:
                    base = _job_without_identity(
                        phase, lane, arm_name, block_count, schedule)
                    jobs.append({**base, "job_id": canonical_sha256(base)})
    if len({job["job_id"] for job in jobs}) != len(jobs):
        raise CampaignError("frozen job identities are not unique")
    order_seed = _SEEDS["order"][phase]["value"].to_bytes(8, "big")
    jobs.sort(key=lambda job: (
        hashlib.sha256(
            order_seed + job["job_id"].encode("ascii")).digest(),
        job["job_id"],
    ))
    return tuple(
        {**job, "order": order} for order, job in enumerate(jobs))


def pre_cpu_job_list_hash(phase, survivor=None):
    return canonical_sha256(list(build_pre_cpu_jobs(phase, survivor)))


def _cell_contract(phase, survivor=None):
    private = _SEEDS["_private"]
    if phase == "training":
        return {
            "phase": phase,
            "arms": [arm["provisional_id"] for arm in ARMS],
            "K_values": list(K_VALUES),
            "schedules": list(SCHEDULES),
            "lanes": [{
                "name": "training",
                "construction_seed_derivation": CONSTRUCTION_DERIVED,
                "paired_roots": list(private["training_roots"]),
                "paired_loss_seeds": list(private["training_losses"]),
            }],
        }
    if phase == "sealed" and survivor in ARM_BY_NAME:
        return {
            "phase": phase,
            "arms": [ARM_BY_NAME[survivor]["provisional_id"]],
            "K_values": list(K_VALUES),
            "schedules": list(SCHEDULES),
            "lanes": [
                {
                    "name": "random",
                    "construction_seed_derivation": CONSTRUCTION_DERIVED,
                    "paired_roots": list(
                        private["sealed_random_roots"]),
                    "paired_loss_seeds": list(
                        private["sealed_random_losses"]),
                },
                {
                    "name": "production",
                    "construction_seed_derivation": CONSTRUCTION_FIXED,
                    "roots_by_K": {
                        str(block_count): _SEEDS["construction"][
                            "sealed_production"]["roots_by_K"][
                                str(block_count)]["root"]
                        for block_count in K_VALUES
                    },
                    "loss_seeds": list(
                        private["sealed_production_losses"]),
                },
            ],
        }
    raise CampaignError("invalid phase/survivor for cell contract")


def frozen_roster():
    return {
        "K_values": list(K_VALUES),
        "schedules": list(SCHEDULES),
        "arms": [dict(arm) for arm in ARMS],
        "repair_policy_sha256": REPAIR_POLICY_SHA256,
        "native_protocol": REQUIRED_REPAIRTIMING_PROTOCOL,
        "native_schema": REQUIRED_REPAIRTIMING_SCHEMA,
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
        "correctness_prelaunch": {
            "block_bytes": list(CORRECTNESS_BLOCK_BYTES),
            "tails": ["full", "partial"],
            "rule": "not-a-campaign-roster-axis",
        },
        "phases": {
            "training": {
                "arms": 2,
                "lanes": ["training"],
                "replicates": TRAINING_REPLICATES,
            },
            "sealed": {
                "arms": 1,
                "lanes": ["random", "production"],
                "replicates_per_lane": SEALED_REPLICATES,
                "winner_source":
                    "authenticated-sole-training-decision-no-CLI-choice",
            },
        },
        "seed_disjointness": seed_disjointness_proof(),
        "selection_policy": selection_policy(),
    }


def frozen_roster_sha256():
    return canonical_sha256(frozen_roster())


def validate_frozen_contract():
    training = build_pre_cpu_jobs("training")
    if len(training) != 1188:
        raise CampaignError("training population is not 1,188 jobs")
    for survivor in ARM_BY_NAME:
        sealed = build_pre_cpu_jobs("sealed", survivor)
        if len(sealed) != 1188:
            raise CampaignError("sealed population is not 1,188 jobs")

    training_cells = len(training) * TRAINING_REPLICATES
    training_selectors = (
        len(ARMS) * len(K_VALUES) * TRAINING_REPLICATES)
    sealed_cells = 1188 * SEALED_REPLICATES
    sealed_selectors = (
        len(K_VALUES) * (SEALED_REPLICATES + 1))
    expected = {
        "training": {
            "jobs": 1188,
            "recovery_cells": 912384,
            "unique_selectors": 152064,
            "attempt_rows": 7299072,
            "native_rows_per_job": 92928,
            "native_rows": 110398464,
        },
        "sealed": {
            "jobs": 1188,
            "recovery_cells": 304128,
            "unique_selectors": 25443,
            "attempt_rows": 2433024,
            "native_rows_per_job": 30976,
            "native_rows": 36799488,
        },
    }
    actual = {
        "training": {
            "jobs": len(training),
            "recovery_cells": training_cells,
            "unique_selectors": training_selectors,
            "attempt_rows": training_cells * REPAIR_ATTEMPTS,
            "native_rows_per_job": TRAINING_REPLICATES * 121,
            "native_rows": training_cells * 121,
        },
        "sealed": {
            "jobs": 1188,
            "recovery_cells": sealed_cells,
            "unique_selectors": sealed_selectors,
            "attempt_rows": sealed_cells * REPAIR_ATTEMPTS,
            "native_rows_per_job": SEALED_REPLICATES * 121,
            "native_rows": sealed_cells * 121,
        },
    }
    if actual != expected:
        raise CampaignError(
            f"frozen campaign cardinalities changed: {actual!r}")
    policy = selection_policy()
    sealed_jobs = {
        survivor: pre_cpu_job_list_hash("sealed", survivor)
        for survivor in ARM_BY_NAME
    }
    sealed_cells_hash = {
        survivor: canonical_sha256(_cell_contract("sealed", survivor))
        for survivor in ARM_BY_NAME
    }
    evidence = {
        "frozen_roster_sha256": frozen_roster_sha256(),
        "policy_sha256": canonical_sha256(policy),
        "seed_disjointness_sha256":
            canonical_sha256(seed_disjointness_proof()),
        "training_job_list_sha256":
            canonical_sha256(list(training)),
        "training_cell_set_sha256":
            canonical_sha256(_cell_contract("training")),
        "sealed_job_list_sha256": sealed_jobs,
        "sealed_cell_set_sha256": sealed_cells_hash,
        "cardinalities": expected,
    }
    evidence["result_free_plan_sha256"] = canonical_sha256({
        "schema": CAMPAIGN_SCHEMA,
        "frozen_roster_sha256": evidence["frozen_roster_sha256"],
        "policy_sha256": evidence["policy_sha256"],
        "seed_disjointness_sha256":
            evidence["seed_disjointness_sha256"],
        "training_job_list_sha256":
            evidence["training_job_list_sha256"],
        "training_cell_set_sha256":
            evidence["training_cell_set_sha256"],
        "sealed_job_list_sha256": sealed_jobs,
        "sealed_cell_set_sha256": sealed_cells_hash,
        "cardinalities": expected,
    })
    return evidence


def make_plan(phase, survivor=None):
    validation = validate_frozen_contract()
    jobs = build_pre_cpu_jobs(phase, survivor)
    return {
        "schema": CAMPAIGN_SCHEMA,
        "phase": phase,
        "survivor": survivor,
        "frozen_roster": frozen_roster(),
        **validation,
        "pre_cpu_job_list_sha256":
            canonical_sha256(list(jobs)),
        "cell_set_sha256":
            canonical_sha256(_cell_contract(phase, survivor)),
        "pre_cpu_jobs": list(jobs),
    }


def make_result_free_plan():
    """Prebind training and both possible sealed templates without choice."""
    validation = validate_frozen_contract()
    return {
        "schema": "wirehair.wh2.rv4a.result-free-plan.v1",
        "campaign_schema": CAMPAIGN_SCHEMA,
        "frozen_roster": frozen_roster(),
        **validation,
        "training": {
            "phase": "training",
            "survivor": None,
            "pre_cpu_job_list_sha256":
                validation["training_job_list_sha256"],
            "cell_set_sha256":
                validation["training_cell_set_sha256"],
            "cardinalities": validation["cardinalities"]["training"],
        },
        "sealed_templates": {
            arm_name: {
                "phase": "sealed",
                "survivor": arm_name,
                "pre_cpu_job_list_sha256":
                    validation["sealed_job_list_sha256"][arm_name],
                "cell_set_sha256":
                    validation["sealed_cell_set_sha256"][arm_name],
                "cardinalities": validation["cardinalities"]["sealed"],
            }
            for arm_name in sorted(ARM_BY_NAME)
        },
        "sealed_selection":
            "authenticated-complete-training-decision-only-no-candidate-"
            "or-lane-override",
    }


def _fsync_directory(directory):
    descriptor = os.open(str(directory), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _commit_no_replace(temporary, path):
    """Atomically publish ``temporary`` without replacing ``path``."""
    try:
        os.link(temporary, path, follow_symlinks=False)
    except FileExistsError:
        raise CampaignError(f"refusing to replace existing file {path}")
    os.unlink(temporary)
    _fsync_directory(path.parent)


def atomic_json(path, value):
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


def atomic_gzip_json(path, value, *, limits=None):
    """Create one deterministic gzip JSON member under optional hard caps."""
    path = Path(path)
    payload = canonical_json_bytes(value)
    if limits is not None:
        if (
            not isinstance(limits, dict) or
            set(limits) != {"compressed", "uncompressed"} or
            any(type(limits.get(name)) is not int or limits[name] < 1
                for name in ("compressed", "uncompressed"))
        ):
            raise CampaignError("gzip publication limits are invalid")
        if len(payload) > limits["uncompressed"]:
            raise CampaignError(
                "result envelope exceeds frozen evidence cap")
    compressed = gzip.compress(payload, compresslevel=6, mtime=0)
    if limits is not None and len(compressed) > limits["compressed"]:
        raise CampaignError("result envelope exceeds frozen evidence cap")
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


def _read_stable_bytes(path, *, byte_limit):
    """Read one regular file through a retained no-follow descriptor."""
    path = Path(path)
    if type(byte_limit) is not int or byte_limit < 0:
        raise CampaignError("stable-read byte limit is invalid")
    flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NONBLOCK", 0)
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    stable = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    descriptor = None
    try:
        descriptor = os.open(str(path), flags)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode) or
            before.st_size > byte_limit
        ):
            raise CampaignError(
                f"evidence is not regular or exceeds byte limit: {path}")
        chunks = []
        count = 0
        while True:
            chunk = os.read(descriptor, min(
                1024 * 1024, byte_limit - count + 1))
            if not chunk:
                break
            count += len(chunk)
            if count > byte_limit:
                raise CampaignError(f"evidence exceeds byte limit: {path}")
            chunks.append(chunk)
        after = os.fstat(descriptor)
        current = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise CampaignError(f"could not stably read {path}: {error}")
    finally:
        if descriptor is not None:
            os.close(descriptor)
    if (
        any(getattr(before, name) != getattr(after, name)
            for name in stable) or
        any(getattr(after, name) != getattr(current, name)
            for name in stable) or
        count != after.st_size
    ):
        raise CampaignError(f"evidence changed while read: {path}")
    return b"".join(chunks)


def read_canonical_json(path, *, byte_limit=64 * 1024 * 1024):
    path = Path(path)
    payload = _read_stable_bytes(path, byte_limit=byte_limit)
    return _strict_json_payload(payload, str(path)), sha256_bytes(payload)


class _EvidenceByteLimitExceeded(Exception):
    pass


class _CappedHashReader:
    def __init__(self, source, limit):
        self.source = source
        self.limit = limit
        self.count = 0
        self.digest = hashlib.sha256()

    def read(self, size=-1):
        remaining_with_probe = self.limit - self.count + 1
        if size is None or size < 0 or size > remaining_with_probe:
            size = remaining_with_probe
        data = self.source.read(size)
        if self.count + len(data) > self.limit:
            raise _EvidenceByteLimitExceeded
        self.count += len(data)
        self.digest.update(data)
        return data


def _result_evidence_byte_limits(job):
    if not isinstance(job, dict):
        raise CampaignError("result evidence job is malformed")
    replicates = job.get("replicates")
    warmups = job.get("warmup_replicates")
    if (
        type(replicates) is not int or not 1 <= replicates <= 768 or
        type(warmups) is not int or warmups < 0
    ):
        raise CampaignError("result evidence replicate count is invalid")
    cells = replicates + warmups
    uncompressed = (
        RESULT_UNCOMPRESSED_FIXED_BYTES
        + cells * RESULT_UNCOMPRESSED_BYTES_PER_CELL)
    # zlib's compressBound formula plus a conservative 64-byte gzip/header
    # allowance.  The compressed cap is derived from, and can never be chosen
    # by, the receipt.  It remains valid even for incompressible canonical
    # evidence rather than assuming an empirical compression ratio.
    compressed = (
        uncompressed + (uncompressed >> 12) +
        (uncompressed >> 14) + (uncompressed >> 25) + 64)
    limits = {
        "compressed": compressed,
        "uncompressed": uncompressed,
    }
    if repair is not None:
        parser_cell_limit = getattr(
            repair,
            "REPAIRTIMING_RECEIPT_CELL_CANONICAL_BYTE_LIMIT",
            None,
        )
        parser_fixed_limit = getattr(
            repair,
            "REPAIRTIMING_RECEIPT_FIXED_CANONICAL_BYTE_LIMIT",
            None,
        )
        if (
            type(parser_cell_limit) is not int or parser_cell_limit < 1 or
            parser_cell_limit > RESULT_UNCOMPRESSED_BYTES_PER_CELL or
            type(parser_fixed_limit) is not int or parser_fixed_limit < 1 or
            parser_fixed_limit > RESULT_UNCOMPRESSED_FIXED_BYTES
        ):
            raise CampaignError(
                "v3 parser receipt cap exceeds campaign cap")
    return limits


def read_canonical_gzip_json(
        path, *, compressed_limit, uncompressed_limit):
    """Read exactly one canonical JSON gzip member under hard byte caps."""
    path = Path(path)
    if (
        type(compressed_limit) is not int or compressed_limit < 1 or
        type(uncompressed_limit) is not int or uncompressed_limit < 1
    ):
        raise CampaignError("gzip evidence byte limits are invalid")
    stable = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    try:
        flags = (
            os.O_RDONLY | os.O_CLOEXEC |
            getattr(os, "O_NONBLOCK", 0)
        )
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(str(path), flags)
        with os.fdopen(descriptor, "rb") as source:
            before = os.fstat(source.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise CampaignError(
                    f"gzip evidence is not a regular file: {path}")
            if before.st_size > compressed_limit:
                raise CampaignError(
                    f"gzip evidence exceeds compressed byte limit: {path}")
            reader = _CappedHashReader(source, compressed_limit)
            chunks = []
            uncompressed_bytes = 0
            archive = zlib.decompressobj(16 + zlib.MAX_WBITS)
            while not archive.eof:
                compressed_chunk = reader.read(64 * 1024)
                if not compressed_chunk:
                    raise CampaignError(
                        f"gzip evidence is truncated: {path}")
                pending = compressed_chunk
                while pending and not archive.eof:
                    remaining_with_probe = (
                        uncompressed_limit - uncompressed_bytes + 1)
                    chunk = archive.decompress(
                        pending,
                        min(64 * 1024, remaining_with_probe),
                    )
                    pending = archive.unconsumed_tail
                    uncompressed_bytes += len(chunk)
                    if uncompressed_bytes > uncompressed_limit:
                        raise CampaignError(
                            "gzip evidence exceeds uncompressed byte "
                            f"limit: {path}")
                    chunks.append(chunk)
                if archive.unused_data:
                    raise CampaignError(
                        f"gzip evidence has trailing data: {path}")
            if reader.read(1):
                raise CampaignError(
                    f"gzip evidence has trailing data: {path}")
            after = os.fstat(source.fileno())
        current = os.stat(path, follow_symlinks=False)
        if (
            any(getattr(before, name) != getattr(after, name)
                for name in stable) or
            any(getattr(after, name) != getattr(current, name)
                for name in stable) or
            reader.count != after.st_size
        ):
            raise CampaignError(
                f"gzip evidence changed while being read: {path}")
        payload = b"".join(chunks)
    except _EvidenceByteLimitExceeded:
        raise CampaignError(
            f"gzip evidence exceeds compressed byte limit: {path}")
    except (OSError, zlib.error) as error:
        raise CampaignError(f"could not read gzip evidence {path}: {error}")
    return (
        _strict_json_payload(payload, str(path)),
        {
            "compressed_sha256": reader.digest.hexdigest(),
            "uncompressed_sha256": sha256_bytes(payload),
            "compressed_bytes": reader.count,
            "uncompressed_bytes": uncompressed_bytes,
        },
    )


def read_canonical_gzip_json_lines(
        path, *, compressed_limit, uncompressed_limit, row_limit):
    """Read one canonical JSONL gzip member and authenticate every row."""
    path = Path(path)
    if (
        type(row_limit) is not int or row_limit < 0 or
        type(compressed_limit) is not int or compressed_limit < 1 or
        type(uncompressed_limit) is not int or uncompressed_limit < 1
    ):
        raise CampaignError("JSONL evidence limits are invalid")
    # Reuse the exact single-member/cap logic by decoding a temporary
    # canonical JSON string would alter evidence hashes, so retain the same
    # bounded zlib loop directly.
    flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NONBLOCK", 0)
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    stable = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    try:
        descriptor = os.open(str(path), flags)
        with os.fdopen(descriptor, "rb") as source:
            before = os.fstat(source.fileno())
            if (
                not stat.S_ISREG(before.st_mode) or
                before.st_size > compressed_limit
            ):
                raise CampaignError(
                    f"JSONL gzip evidence is invalid or oversized: {path}")
            reader = _CappedHashReader(source, compressed_limit)
            archive = zlib.decompressobj(16 + zlib.MAX_WBITS)
            chunks = []
            total = 0
            while not archive.eof:
                compressed = reader.read(64 * 1024)
                if not compressed:
                    raise CampaignError(
                        f"JSONL gzip evidence is truncated: {path}")
                pending = compressed
                while pending and not archive.eof:
                    remaining = uncompressed_limit - total + 1
                    chunk = archive.decompress(
                        pending, min(64 * 1024, remaining))
                    pending = archive.unconsumed_tail
                    total += len(chunk)
                    if total > uncompressed_limit:
                        raise CampaignError(
                            "JSONL gzip evidence exceeds uncompressed cap")
                    chunks.append(chunk)
                if archive.unused_data:
                    raise CampaignError(
                        f"JSONL gzip evidence has trailing data: {path}")
            if reader.read(1):
                raise CampaignError(
                    f"JSONL gzip evidence has trailing data: {path}")
            after = os.fstat(source.fileno())
        current = os.stat(path, follow_symlinks=False)
    except _EvidenceByteLimitExceeded:
        raise CampaignError(
            f"JSONL gzip evidence exceeds compressed cap: {path}")
    except (OSError, zlib.error) as error:
        raise CampaignError(
            f"could not read JSONL gzip evidence {path}: {error}")
    if (
        any(getattr(before, name) != getattr(after, name)
            for name in stable) or
        any(getattr(after, name) != getattr(current, name)
            for name in stable) or
        reader.count != after.st_size
    ):
        raise CampaignError(f"JSONL gzip evidence changed while read: {path}")
    payload = b"".join(chunks)
    if payload and not payload.endswith(b"\n"):
        raise CampaignError("JSONL evidence lacks final LF")
    lines = payload.splitlines(keepends=True)
    if len(lines) > row_limit:
        raise CampaignError("JSONL evidence exceeds row cap")
    rows = [
        _strict_json_payload(line, f"{path}:row:{index}")
        for index, line in enumerate(lines)
    ]
    return rows, {
        "compressed_sha256": reader.digest.hexdigest(),
        "uncompressed_sha256": sha256_bytes(payload),
        "compressed_bytes": reader.count,
        "uncompressed_bytes": len(payload),
        "rows": len(rows),
    }


def _stable_file_binding(
        path, *, executable=False, byte_limit=PYTHON_SOURCE_BYTE_LIMIT):
    if type(byte_limit) is not int or byte_limit < 1:
        raise CampaignError("runtime file byte limit is invalid")
    path = Path(path).resolve()
    flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NONBLOCK", 0)
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(str(path), flags)
        try:
            before = os.fstat(descriptor)
            if (
                not stat.S_ISREG(before.st_mode) or
                before.st_size > byte_limit
            ):
                raise CampaignError(
                    f"runtime file is nonregular or exceeds byte limit: {path}")
            digest = hashlib.sha256()
            size = 0
            while True:
                chunk = os.read(descriptor, 1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
                size += len(chunk)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        current = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise CampaignError(f"could not bind runtime file {path}: {error}")
    stable = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    if (
        any(getattr(before, name) != getattr(after, name)
            for name in stable) or
        any(getattr(after, name) != getattr(current, name)
            for name in stable) or
        not stat.S_ISREG(after.st_mode) or size != after.st_size or
        (executable and after.st_mode & 0o111 == 0)
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
        "byte_limit": byte_limit,
    }


def _thermal_binding(path):
    path = Path(path).resolve()
    try:
        descriptor = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise CampaignError(f"could not bind thermal source {path}: {error}")
    if not stat.S_ISREG(descriptor.st_mode):
        raise CampaignError(f"thermal source is not regular: {path}")
    return {
        "path": str(path),
        "device": descriptor.st_dev,
        "inode": descriptor.st_ino,
        "mode": descriptor.st_mode,
    }


def _require_v3_parser():
    if repair is None:
        raise CampaignError("repair_timing_codec v3 parser is unavailable")
    if (
        getattr(repair, "REPAIRTIMING_PROTOCOL_V3", None)
            != REQUIRED_REPAIRTIMING_PROTOCOL or
        getattr(repair, "REPAIRTIMING_SCHEMA_V3", None)
            != REQUIRED_REPAIRTIMING_SCHEMA or
        getattr(repair, "REPAIRTIMING_ROWS_PER_CELL", None) != 121
    ):
        raise CampaignError("repair_timing_codec protocol is not frozen v3")
    parser_path = Path(repair.__file__).resolve()
    context_path = Path(peel_codec.__file__).resolve()
    if (
        parser_path != (TOOLS / "repair_timing_codec.py").resolve() or
        context_path != (TOOLS / "peel_codec.py").resolve()
    ):
        raise CampaignError(
            "local v3 consumers were not loaded from repository source")
    return repair


def runtime_bindings(bench_path, thermal_source):
    parser_module = _require_v3_parser()
    resolved_bench = Path(bench_path).resolve()
    build_configuration = resolved_bench.parent.parent / "CMakeCache.txt"
    sources = {
        "native_benchmark": (
            REPOSITORY / "codec/WirehairV2Bench.cpp",
            NATIVE_SOURCE_BYTE_LIMIT),
        "native_repairtiming": (
            REPOSITORY / "codec/WirehairV2RepairTiming.inc",
            NATIVE_SOURCE_BYTE_LIMIT),
        "repair_contract": (
            REPOSITORY / "codec/WirehairV2Repair.cpp",
            NATIVE_SOURCE_BYTE_LIMIT),
        "repair_contract_header": (
            REPOSITORY / "codec/WirehairV2Repair.h",
            NATIVE_SOURCE_BYTE_LIMIT),
        "native_cli_test": (
            REPOSITORY / "codec/V2BenchCliTest.cmake",
            NATIVE_SOURCE_BYTE_LIMIT),
        "native_repairtiming_cli_test": (
            REPOSITORY / "codec/V2RepairTimingCliTest.cmake",
            NATIVE_SOURCE_BYTE_LIMIT),
        "build_policy_e2e": (
            REPOSITORY / "test/cmake/RunBuildPolicyE2E.cmake",
            BUILD_SOURCE_BYTE_LIMIT),
        "codec_build": (
            REPOSITORY / "codec/CMakeLists.txt",
            BUILD_SOURCE_BYTE_LIMIT),
        "root_build": (
            REPOSITORY / "CMakeLists.txt",
            BUILD_SOURCE_BYTE_LIMIT),
        "parser": (
            Path(parser_module.__file__), PYTHON_SOURCE_BYTE_LIMIT),
        "context_tool": (
            Path(peel_codec.__file__), PYTHON_SOURCE_BYTE_LIMIT),
        "runner": (Path(__file__), PYTHON_SOURCE_BYTE_LIMIT),
        "parser_test": (
            REPOSITORY / "tools/test_repair_timing_tools.py",
            PYTHON_SOURCE_BYTE_LIMIT),
        "runner_test": (
            REPOSITORY / "bench/test_wh2_rv4a_campaign.py",
            PYTHON_SOURCE_BYTE_LIMIT),
        "parallel_verifier": (
            REPOSITORY / "bench/wh2_rv4a_parallel_verify.py",
            PYTHON_SOURCE_BYTE_LIMIT),
        "parallel_verifier_test": (
            REPOSITORY / "bench/test_wh2_rv4a_parallel_verify.py",
            PYTHON_SOURCE_BYTE_LIMIT),
        "build_configuration": (
            build_configuration, BUILD_CONFIGURATION_BYTE_LIMIT),
    }
    bindings = {
        "benchmark": _stable_file_binding(
            resolved_bench,
            executable=True,
            byte_limit=BENCHMARK_BINARY_BYTE_LIMIT,
        ),
        "sources": {
            name: _stable_file_binding(path, byte_limit=limit)
            for name, (path, limit) in sources.items()
        },
        "thermal": _thermal_binding(thermal_source),
    }
    _validate_loaded_source_hashes(bindings, parser_module)
    return bindings


def _validate_loaded_source_hashes(bindings, parser_module):
    """Reject a file mutation between source compilation and binding."""
    loaded_hashes = {
        "parser": getattr(
            parser_module, "__rv4a_source_sha256__", None),
        "context_tool": getattr(
            peel_codec, "__rv4a_source_sha256__", None),
    }
    if _BOOTSTRAP_RUNNER_SOURCE_SHA256 is not None:
        loaded_hashes["runner"] = _BOOTSTRAP_RUNNER_SOURCE_SHA256
    if any(
            not _is_sha256(loaded_hashes[name]) or
            bindings["sources"][name]["sha256"] != loaded_hashes[name]
            for name in loaded_hashes):
        raise CampaignError(
            "runtime source changed after source-forced bootstrap")
    return bindings


def _require_source_forced_outcome_runner(expected_sha256=None):
    """Fail closed unless outcome execution came through the exact-source CLI."""
    if not _is_sha256(_BOOTSTRAP_RUNNER_SOURCE_SHA256):
        raise CampaignError(
            "outcome execution requires the source-forced campaign runner")
    if (
        expected_sha256 is not None and
        (
            not _is_sha256(expected_sha256) or
            _BOOTSTRAP_RUNNER_SOURCE_SHA256 != expected_sha256
        )
    ):
        raise CampaignError(
            "executing campaign runner does not match outcome manifest")


def _validate_runtime_bindings_schema(bindings):
    file_fields = {
        "path", "device", "inode", "mode", "size",
        "mtime_ns", "ctime_ns", "sha256", "byte_limit",
    }
    source_names = {
        "native_benchmark", "native_repairtiming", "repair_contract",
        "repair_contract_header", "native_cli_test",
        "native_repairtiming_cli_test", "build_policy_e2e",
        "codec_build", "root_build", "parser", "context_tool", "runner",
        "parser_test", "runner_test", "parallel_verifier",
        "parallel_verifier_test", "build_configuration",
    }
    if (
        not isinstance(bindings, dict) or
        set(bindings) != {"benchmark", "sources", "thermal"} or
        not isinstance(bindings["sources"], dict) or
        set(bindings["sources"]) != source_names
    ):
        raise CampaignError("runtime binding schema is invalid")
    for item in (bindings["benchmark"], *bindings["sources"].values()):
        if (
            not isinstance(item, dict) or set(item) != file_fields or
            not isinstance(item.get("path"), str) or
            not Path(item["path"]).is_absolute() or
            any(type(item.get(field)) is not int or item[field] < 0
                for field in (
                    "device", "inode", "mode", "size",
                    "mtime_ns", "ctime_ns", "byte_limit")) or
            item["size"] > item["byte_limit"] or
            not _is_sha256(item.get("sha256"))
        ):
            raise CampaignError("runtime file binding is invalid")
    thermal = bindings["thermal"]
    if (
        not isinstance(thermal, dict) or
        set(thermal) != {"path", "device", "inode", "mode"} or
        not isinstance(thermal["path"], str) or
        not Path(thermal["path"]).is_absolute() or
        any(type(thermal.get(field)) is not int or thermal[field] < 0
            for field in ("device", "inode", "mode"))
    ):
        raise CampaignError("thermal binding is invalid")
    return bindings


_RUNTIME_PIN_SOURCE_CLASSES = {
    "native_benchmark": "native",
    "native_repairtiming": "native",
    "repair_contract": "native",
    "repair_contract_header": "native",
    "native_cli_test": "native",
    "native_repairtiming_cli_test": "native",
    "build_policy_e2e": "build",
    "codec_build": "build",
    "root_build": "build",
    "parser": "python",
    "context_tool": "python",
    "runner": "python",
    "parser_test": "python",
    "runner_test": "python",
    "parallel_verifier": "python",
    "parallel_verifier_test": "python",
    "build_configuration": "configuration",
}


def _runtime_pin_projection(bindings):
    _validate_runtime_bindings_schema(bindings)
    files = {"benchmark": bindings["benchmark"], **bindings["sources"]}
    return {
        name: {
            "class": (
                "binary" if name == "benchmark"
                else _RUNTIME_PIN_SOURCE_CLASSES[name]
            ),
            "sha256": item["sha256"],
            "size": item["size"],
            "byte_limit": item["byte_limit"],
        }
        for name, item in sorted(files.items())
    }


def make_runtime_pin_record(bindings):
    """Return the exact pre-outcome source/binary pin for these bindings."""
    validation = validate_frozen_contract()
    return {
        "schema": RUNTIME_PIN_SCHEMA,
        "campaign_schema": CAMPAIGN_SCHEMA,
        "repairtiming_protocol": REQUIRED_REPAIRTIMING_PROTOCOL,
        "repairtiming_schema": REQUIRED_REPAIRTIMING_SCHEMA,
        "frozen_roster_sha256": validation["frozen_roster_sha256"],
        "policy_sha256": validation["policy_sha256"],
        "result_free_plan_sha256":
            validation["result_free_plan_sha256"],
        "runtime_files": _runtime_pin_projection(bindings),
        "thermal": dict(bindings["thermal"]),
    }


def make_preflight_pin(bench_path, thermal_source):
    """Hash the complete runtime before any campaign directory or outcome."""
    bindings = runtime_bindings(bench_path, thermal_source)
    return make_runtime_pin_record(bindings)


def _validate_preflight_pin(path, expected_sha256, bindings):
    if not _is_sha256(expected_sha256):
        raise CampaignError("trusted preflight pin digest is invalid")
    pin, digest = read_canonical_json(
        path, byte_limit=RUNTIME_PIN_BYTE_LIMIT)
    expected = make_runtime_pin_record(bindings)
    if (
        digest != expected_sha256 or
        canonical_json_bytes(pin) != canonical_json_bytes(expected)
    ):
        raise CampaignError(
            "preflight pin does not match the complete runtime")
    return digest


def check_runtime_bindings(expected, *, full_hash):
    _validate_runtime_bindings_schema(expected)
    if _thermal_binding(expected["thermal"]["path"]) != expected["thermal"]:
        raise CampaignError("thermal source identity changed")
    named = {"benchmark": expected["benchmark"], **expected["sources"]}
    for name, item in named.items():
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
                path,
                executable=(name == "benchmark"),
                byte_limit=item["byte_limit"],
            )
            if rebound != item:
                raise CampaignError(f"runtime content changed for {name}")


def _require_normal_priority():
    if not hasattr(os, "getpriority"):
        raise CampaignError("normal-priority validation is unavailable")
    priority = os.getpriority(os.PRIO_PROCESS, 0)
    if priority != 0:
        raise CampaignError(
            f"campaign workers must run at nice 0, got {priority}")


def _create_fresh_directory(directory):
    directory = Path(directory)
    if directory.exists():
        raise CampaignError(f"refusing existing directory {directory}")
    try:
        directory.mkdir(parents=True)
    except FileExistsError:
        raise CampaignError(f"refusing existing directory {directory}")
    _fsync_directory(directory.parent)
    return directory


def expected_request(job):
    return {
        "block_count": job["K"],
        "block_bytes": job["block_bytes"],
        "repair_arm": job["arm_native_name"],
        "repair_policy": "repair-v1",
        "dispatch_profile": "dispatch-v1",
        "construction_seed": job["construction_seed_base"],
        "construction_seed_derivation":
            job["construction_seed_derivation"],
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
    lane = job.get("lane")
    arm = job.get("arm")
    block_count = job.get("K")
    schedule = job.get("schedule")
    expected_lanes = (
        ("training",) if phase == "training"
        else ("random", "production") if phase == "sealed"
        else ()
    )
    if (
        lane not in expected_lanes or arm not in ARM_BY_NAME or
        block_count not in K_VALUES or schedule not in SCHEDULES
    ):
        raise CampaignError("job is outside the frozen roster")
    expected = _job_without_identity(
        phase, lane, arm, block_count, schedule)
    order = job.get("order")
    if (
        type(order) is not int or order < 0 or
        not peel_codec._same_typed_json(job, {
            **expected,
            "job_id": canonical_sha256(expected),
            "order": order,
        })
    ):
        raise CampaignError("job differs from its frozen request")
    request = expected_request(job)
    if repair is not None:
        parser_module = _require_v3_parser()
        parser_module.validate_repairtiming_dimensions(**request)
    return request


def _expected_assignment(index, job, cpus):
    return {
        "order": index,
        "job_id": job["job_id"],
        "cpu": cpus[index % len(cpus)],
        "job_file": f"jobs/{index:04d}.json",
        "output": f"results/{index:04d}.json.gz",
        "log": f"logs/{index:04d}.log",
    }


def _build_manifest(
        phase, survivor, jobs, cpus, bindings, training_contract,
        preflight_pin_sha256):
    plan = make_plan(phase, survivor)
    return {
        **plan,
        "created_unix_ns": time.time_ns(),
        "repository": str(REPOSITORY),
        "workers": len(cpus),
        "worker_cpus": cpus,
        "priority": "normal-nice-0",
        "preflight_pin_sha256": preflight_pin_sha256,
        "runtime_bindings": bindings,
        "training_contract": training_contract,
        "assignments": [
            _expected_assignment(index, job, cpus)
            for index, job in enumerate(jobs)
        ],
    }


def _validate_manifest(manifest):
    runtime_fields = {
        "created_unix_ns", "repository", "workers", "worker_cpus",
        "priority", "preflight_pin_sha256", "runtime_bindings",
        "training_contract", "assignments",
    }
    if (
        not isinstance(manifest, dict) or
        manifest.get("schema") != CAMPAIGN_SCHEMA or
        manifest.get("phase") not in ("training", "sealed")
    ):
        raise CampaignError("campaign manifest schema is invalid")
    phase = manifest["phase"]
    survivor = manifest.get("survivor")
    if (
        (phase == "training" and survivor is not None) or
        (phase == "sealed" and survivor not in ARM_BY_NAME)
    ):
        raise CampaignError("campaign manifest survivor is invalid")
    expected_plan = make_plan(phase, survivor)
    cpus = manifest.get("worker_cpus")
    if (
        set(manifest) != set(expected_plan) | runtime_fields or
        any(not peel_codec._same_typed_json(
                manifest.get(key), value)
            for key, value in expected_plan.items()) or
        type(manifest.get("created_unix_ns")) is not int or
        manifest["created_unix_ns"] <= 0 or
        manifest.get("repository") != str(REPOSITORY) or
        type(manifest.get("workers")) is not int or
        not 1 <= manifest["workers"] <= MAX_WORKERS or
        not isinstance(cpus, list) or len(cpus) != manifest["workers"] or
        cpus != sorted(set(cpus)) or
        any(type(cpu) is not int or cpu < 0 for cpu in cpus) or
        manifest.get("priority") != "normal-nice-0" or
        not _is_sha256(manifest.get("preflight_pin_sha256")) or
        manifest["preflight_pin_sha256"] != canonical_sha256(
            make_runtime_pin_record(manifest["runtime_bindings"]))
    ):
        raise CampaignError("campaign manifest runtime contract is invalid")
    _validate_runtime_bindings_schema(manifest["runtime_bindings"])
    training_contract = manifest.get("training_contract")
    if phase == "training":
        if training_contract is not None:
            raise CampaignError("training manifest fabricated a parent")
    else:
        _validate_training_contract_shape(
            training_contract, expected_survivor=survivor)
    jobs = build_pre_cpu_jobs(phase, survivor)
    assignments = manifest.get("assignments")
    if (
        not isinstance(assignments, list) or
        len(assignments) != len(jobs)
    ):
        raise CampaignError("campaign assignments are incomplete")
    for index, (assignment, job) in enumerate(zip(assignments, jobs)):
        if not peel_codec._same_typed_json(
                assignment, _expected_assignment(index, job, cpus)):
            raise CampaignError("campaign CPU assignment changed")
    return jobs


def _write_job_files(directory, manifest, manifest_sha256):
    evidence = []
    for assignment, job in zip(
            manifest["assignments"], manifest["pre_cpu_jobs"]):
        record = {
            "schema": JOB_SCHEMA,
            "manifest_sha256": manifest_sha256,
            "runtime_bindings": manifest["runtime_bindings"],
            "job": job,
            "assignment": assignment,
        }
        path = directory / assignment["job_file"]
        digest = atomic_json(path, record)
        evidence.append({
            "path": assignment["job_file"],
            "sha256": digest,
            "bytes": len(canonical_json_bytes(record)),
        })
    return evidence


def _verify_job_files(directory, manifest, manifest_sha256):
    evidence = []
    for assignment, job in zip(
            manifest["assignments"], manifest["pre_cpu_jobs"]):
        path = directory / assignment["job_file"]
        record, digest = read_canonical_json(path, byte_limit=4 * 1024 * 1024)
        expected = {
            "schema": JOB_SCHEMA,
            "manifest_sha256": manifest_sha256,
            "runtime_bindings": manifest["runtime_bindings"],
            "job": job,
            "assignment": assignment,
        }
        if not peel_codec._same_typed_json(record, expected):
            raise CampaignError(f"job file changed: {path}")
        evidence.append({
            "path": assignment["job_file"],
            "sha256": digest,
            "bytes": len(canonical_json_bytes(record)),
        })
    return evidence


def _validate_worker_record(record):
    if (
        not isinstance(record, dict) or
        set(record) != {
            "schema", "manifest_sha256", "runtime_bindings",
            "job", "assignment",
        } or
        record.get("schema") != JOB_SCHEMA or
        not _is_sha256(record.get("manifest_sha256"))
    ):
        raise CampaignError("worker job record is malformed")
    _validate_runtime_bindings_schema(record["runtime_bindings"])
    job = record["job"]
    validate_job(job)
    assignment = record["assignment"]
    if (
        not isinstance(assignment, dict) or
        set(assignment) != {
            "order", "job_id", "cpu", "job_file", "output", "log"
        } or
        type(assignment.get("order")) is not int or
        assignment["order"] != job["order"] or
        assignment.get("job_id") != job["job_id"] or
        type(assignment.get("cpu")) is not int or
        assignment["cpu"] < 0 or
        any(type(assignment.get(name)) is not str
            for name in ("job_file", "output", "log"))
    ):
        raise CampaignError("worker assignment is malformed")
    return job, assignment


def _measurement_as_dict(measurement):
    try:
        receipt = measurement.as_dict()
    except AttributeError:
        raise CampaignError("repairtiming parser returned no receipt")
    if not isinstance(receipt, dict):
        raise CampaignError("repairtiming receipt is not a dictionary")
    return receipt


def run_worker(job_file, output_path):
    """Execute and immediately replay one pinned native v3 job."""
    _require_source_forced_outcome_runner()
    _require_normal_priority()
    parser_module = _require_v3_parser()
    record, unused_digest = read_canonical_json(
        job_file, byte_limit=4 * 1024 * 1024)
    job, assignment = _validate_worker_record(record)
    job_file = Path(job_file).resolve()
    output_path = Path(output_path).resolve()
    directory = job_file.parent.parent
    if (
        job_file != (directory / assignment["job_file"]).resolve() or
        output_path != (directory / assignment["output"]).resolve()
    ):
        raise CampaignError("worker paths do not match assignment")
    try:
        os.sched_setaffinity(0, {assignment["cpu"]})
    except (AttributeError, OSError) as error:
        raise CampaignError(f"could not pin worker CPU: {error}")
    if sorted(os.sched_getaffinity(0)) != [assignment["cpu"]]:
        raise CampaignError("worker CPU affinity does not match assignment")

    bindings = record["runtime_bindings"]
    before = runtime_bindings(
        bindings["benchmark"]["path"], bindings["thermal"]["path"])
    if not peel_codec._same_typed_json(before, bindings):
        raise CampaignError("worker runtime binding changed before launch")
    request = expected_request(job)
    context = peel_codec.make_paired_context_config(
        bindings["thermal"]["path"],
        job["thermal_sampling_interval_ms"],
        cache_state=request["cache_state"],
        evict_bytes=request["evict_bytes"],
    )
    started = time.time_ns()
    measurement = parser_module.repairtiming_probe(
        bindings["benchmark"]["path"],
        context=context,
        **request,
    )
    finished = time.time_ns()
    if getattr(measurement, "protocol", None) != REQUIRED_REPAIRTIMING_PROTOCOL:
        raise CampaignError("native result is not repairtiming v3")
    receipt = _measurement_as_dict(measurement)
    replayed = parser_module.replay_repairtiming_receipt(
        receipt, expected_request=request)
    if (
        getattr(replayed, "protocol", None) != REQUIRED_REPAIRTIMING_PROTOCOL or
        not peel_codec._same_typed_json(
            _measurement_as_dict(replayed), receipt)
    ):
        raise CampaignError("immediate parser replay changed the receipt")
    context_receipt = receipt.get("context")
    bound = (
        context_receipt.get("bound")
        if isinstance(context_receipt, dict) else None
    )
    thermal = (
        context_receipt.get("thermal")
        if isinstance(context_receipt, dict) else None
    )
    if (
        not isinstance(bound, dict) or
        bound.get("cpu_affinity") != [assignment["cpu"]] or
        bound.get("thermal_device") != bindings["thermal"]["device"] or
        bound.get("thermal_inode") != bindings["thermal"]["inode"] or
        not isinstance(thermal, dict) or
        thermal.get("edac_ce_max") != 0 or
        thermal.get("edac_ue_max") != 0 or
        type(thermal.get("valid_rows")) is not int or
        thermal["valid_rows"] < 2
    ):
        raise CampaignError(
            "receipt lacks assigned CPU/thermal/EDAC coverage")
    after = runtime_bindings(
        bindings["benchmark"]["path"], bindings["thermal"]["path"])
    if not peel_codec._same_typed_json(after, before):
        raise CampaignError("worker runtime binding changed during probe")
    envelope = {
        "schema": RESULT_SCHEMA,
        "manifest_sha256": record["manifest_sha256"],
        "job": job,
        "assignment": assignment,
        "wall_started_unix_ns": started,
        "wall_finished_unix_ns": finished,
        "runtime_bindings_before": before,
        "runtime_bindings_after": after,
        "receipt": receipt,
    }
    atomic_gzip_json(
        output_path,
        envelope,
        limits=_result_evidence_byte_limits(job),
    )


def _parent_signal_error(signum):
    return CampaignError(
        "campaign parent received "
        f"{signal.Signals(signum).name}")


def _raise_pending_parent_signal(state):
    for item in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        signal.signal(item, signal.SIG_IGN)
    raise _parent_signal_error(state["pending"])


def _parent_signal_safe_point():
    state = _PARENT_SIGNAL_STATE
    if state is not None and state["pending"] is not None:
        _raise_pending_parent_signal(state)


@contextmanager
def _parent_termination_handlers():
    global _PARENT_SIGNAL_STATE
    if _PARENT_SIGNAL_STATE is not None:
        yield
        return
    state = {"pending": None}
    previous = {}

    def handler(signum, unused_frame):
        if state["pending"] is None:
            state["pending"] = signum

    for item in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        previous[item] = signal.signal(item, handler)
    _PARENT_SIGNAL_STATE = state
    try:
        yield
    finally:
        _PARENT_SIGNAL_STATE = None
        for item, prior in previous.items():
            signal.signal(item, prior)


def _process_group_exists(group):
    try:
        os.killpg(group, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def _pin_worker_before_exec(cpu):
    """Restrict a rolling worker before Python imports execute after exec()."""
    if type(cpu) is not int or cpu < 0:
        raise OSError("rolling worker CPU is invalid")
    os.sched_setaffinity(0, {cpu})
    if sorted(os.sched_getaffinity(0)) != [cpu]:
        raise OSError("rolling worker pre-exec affinity did not take effect")


def _kill_process_groups_with_grace(processes, grace_seconds):
    """Terminate owned worker groups, reap leaders, and prove no survivor."""
    if (
        type(grace_seconds) not in (int, float) or
        isinstance(grace_seconds, bool) or
        not math.isfinite(float(grace_seconds)) or grace_seconds < 0
    ):
        raise CampaignError("worker termination grace is invalid")
    unique = {}
    for process in processes:
        pid = getattr(process, "pid", None)
        if type(pid) is not int or pid <= 0:
            continue
        prior = unique.get(pid)
        if prior is not None and prior is not process:
            raise CampaignError("worker cleanup PID ownership is ambiguous")
        unique[pid] = process
    groups = {
        pid for pid in unique if _process_group_exists(pid)
    }
    for group in sorted(groups):
        try:
            os.killpg(group, signal.SIGTERM)
        except ProcessLookupError:
            pass
    for pid, process in unique.items():
        if pid not in groups and process.poll() is None:
            process.terminate()

    def pending(tracked_groups):
        returncodes = {
            pid: process.poll()
            for pid, process in unique.items()
        }
        live_groups = {
            group for group in tracked_groups
            if _process_group_exists(group)
        }
        live_leaders = {
            pid for pid, returncode in returncodes.items()
            if returncode is None
        }
        return live_groups, live_leaders

    deadline = time.monotonic() + grace_seconds
    live_groups, live_leaders = pending(groups)
    while (
        (live_groups or live_leaders) and
        time.monotonic() < deadline
    ):
        if live_groups or live_leaders:
            time.sleep(0.02)
        live_groups, live_leaders = pending(live_groups)
    for group in sorted(live_groups):
        try:
            os.killpg(group, signal.SIGKILL)
        except ProcessLookupError:
            pass
    for pid in sorted(live_leaders):
        process = unique[pid]
        if pid not in live_groups:
            process.kill()

    deadline = time.monotonic() + WORKER_KILL_GRACE_SECONDS
    live_groups, live_leaders = pending(live_groups)
    while (
        (live_groups or live_leaders) and
        time.monotonic() < deadline
    ):
        time.sleep(0.02)
        live_groups, live_leaders = pending(live_groups)
    for process in unique.values():
        try:
            process.wait(timeout=0)
        except subprocess.TimeoutExpired:
            pass
    live_groups, live_leaders = pending(live_groups)
    if live_groups or live_leaders:
        raise CampaignError(
            "rolling worker cleanup left surviving process groups/leaders: "
            f"groups={sorted(live_groups)!r} "
            f"leaders={sorted(live_leaders)!r}")


def _run_rolling_workers(directory, manifest):
    """Keep one worker active on every CPU until that CPU's queue drains."""
    cpus = manifest["worker_cpus"]
    queues = {
        cpu: [
            assignment for assignment in manifest["assignments"]
            if assignment["cpu"] == cpu
        ]
        for cpu in cpus
    }
    cursors = {cpu: 0 for cpu in cpus}
    active = {}
    owned_processes = {}
    queue_receipts = {
        cpu: {
            "cpu": cpu,
            "orders": [
                assignment["order"] for assignment in queues[cpu]
            ],
            "launches": [],
        }
        for cpu in cpus
    }
    check_count = 0
    checks_running = 0
    pending_idle_checks = 0
    occupancy_histogram = {}
    drained_with_global_unstarted_checks = 0
    completed_groups_confirmed_gone = 0
    first_check = last_check = 0

    def launch(cpu):
        if cpu in active or cursors[cpu] >= len(queues[cpu]):
            raise CampaignError("rolling worker launch state is invalid")
        assignment = queues[cpu][cursors[cpu]]
        cursors[cpu] += 1
        log = open(directory / assignment["log"], "xb")
        try:
            environment = os.environ.copy()
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            command = [
                sys.executable,
                "-B",
                str(Path(__file__).resolve()),
                "_worker",
                "--job-file", str(directory / assignment["job_file"]),
                "--output", str(directory / assignment["output"]),
            ]
            process = subprocess.Popen(
                command,
                stdin=subprocess.DEVNULL,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=environment,
                start_new_session=True,
                preexec_fn=functools.partial(
                    _pin_worker_before_exec, cpu),
            )
        except BaseException:
            log.close()
            raise
        started = time.monotonic_ns()
        receipt = {
            "order": assignment["order"],
            "started_monotonic_ns": started,
            "finished_monotonic_ns": 0,
            "returncode": None,
        }
        queue_receipts[cpu]["launches"].append(receipt)
        active[cpu] = {
            "assignment": assignment,
            "process": process,
            "log": log,
            "receipt": receipt,
        }
        if process.pid is None or process.pid in owned_processes:
            try:
                process.kill()
            finally:
                raise CampaignError(
                    "rolling worker returned an invalid or reused PID")
        owned_processes[process.pid] = process

    def sample():
        nonlocal check_count, checks_running, pending_idle_checks
        nonlocal drained_with_global_unstarted_checks
        nonlocal first_check, last_check
        idle_with_pending = [
            cpu for cpu in cpus
            if cursors[cpu] < len(queues[cpu]) and cpu not in active
        ]
        if idle_with_pending:
            pending_idle_checks += 1
            raise CampaignError(
                "rolling scheduler left assigned CPU work idle: "
                f"{idle_with_pending!r}")
        global_unstarted = sum(
            len(queues[cpu]) - cursors[cpu] for cpu in cpus)
        if global_unstarted and any(
                cpu not in active and cursors[cpu] == len(queues[cpu])
                for cpu in cpus):
            # CPU affinity and the assignment are frozen in the manifest.
            # Consequently a CPU may become idle only after its own immutable
            # queue drains, even if another CPU still has assigned jobs.
            drained_with_global_unstarted_checks += 1
        now = time.monotonic_ns()
        check_runtime_bindings(
            manifest["runtime_bindings"], full_hash=False)
        if first_check == 0:
            first_check = now
        last_check = now
        check_count += 1
        active_count = len(active)
        if active_count:
            checks_running += 1
        key = str(active_count)
        occupancy_histogram[key] = occupancy_histogram.get(key, 0) + 1
        _parent_signal_safe_point()

    try:
        for cpu in cpus:
            _parent_signal_safe_point()
            launch(cpu)
        if len(active) != manifest["workers"]:
            raise CampaignError(
                "rolling scheduler did not fill every startup slot")
        sample()
        while active:
            descendant_groups = []
            finished = []
            for cpu in sorted(active):
                state = active[cpu]
                returncode = state["process"].poll()
                if returncode is None:
                    continue
                state["receipt"]["finished_monotonic_ns"] = \
                    time.monotonic_ns()
                state["receipt"]["returncode"] = returncode
                state["log"].close()
                process = state["process"]
                group_exists = _process_group_exists(process.pid)
                if group_exists:
                    descendant_groups.append((
                        state["assignment"]["order"], process.pid))
                else:
                    owned_processes.pop(process.pid, None)
                    completed_groups_confirmed_gone += 1
                finished.append((cpu, returncode))
            for cpu, unused_returncode in finished:
                del active[cpu]
            failures = [
                (queue_receipts[cpu]["launches"][-1]["order"], returncode)
                for cpu, returncode in finished
                if returncode != 0
            ]
            if failures:
                raise CampaignError(
                    f"rolling campaign worker failed: {failures!r}")
            if descendant_groups:
                raise CampaignError(
                    "rolling campaign worker left descendants: "
                    f"{descendant_groups!r}")
            for cpu, unused_returncode in finished:
                if cursors[cpu] < len(queues[cpu]):
                    launch(cpu)
            sample()
            if active:
                time.sleep(WORKER_POLL_SECONDS)
    except BaseException:
        _kill_process_groups_with_grace(
            list(owned_processes.values()),
            WORKER_TERMINATE_GRACE_SECONDS,
        )
        raise
    finally:
        for state in active.values():
            state["log"].close()
    if owned_processes:
        raise CampaignError("rolling scheduler lost owned worker groups")
    return {
        "scheduler": "per-cpu-rolling-v1",
        "checks": check_count,
        "checks_while_running": checks_running,
        "first_check_monotonic_ns": first_check,
        "last_check_monotonic_ns": last_check,
        "runtime_bindings_sha256":
            canonical_sha256(manifest["runtime_bindings"]),
        "startup_slots": manifest["workers"],
        "pending_idle_checks": pending_idle_checks,
        "drain_policy":
            "idle-only-after-that-cpu-fixed-queue-exhausted-v1",
        "drained_slots_while_global_unstarted_checks":
            drained_with_global_unstarted_checks,
        "completed_groups_confirmed_gone":
            completed_groups_confirmed_gone,
        "occupancy_histogram": occupancy_histogram,
        "cpu_queues": [
            queue_receipts[cpu] for cpu in cpus
        ],
    }


def _verify_empty_logs(directory, manifest):
    evidence = []
    for assignment in manifest["assignments"]:
        path = directory / assignment["log"]
        payload = _read_stable_bytes(path, byte_limit=0)
        if payload:
            raise CampaignError(f"worker log is not empty: {path}")
        evidence.append({
            "path": assignment["log"],
            "sha256": sha256_bytes(payload),
            "bytes": 0,
        })
    return evidence


def _expected_root(job, replicate):
    if job["construction_seed_derivation"] == CONSTRUCTION_DERIVED:
        return construction_seed(job["construction_seed_base"], replicate)
    if job["construction_seed_derivation"] == CONSTRUCTION_FIXED:
        return job["construction_seed_base"]
    raise CampaignError("job has unknown construction derivation")


def _expected_loss(job, replicate):
    if job["loss_seed_derivation"] != LOSS_DERIVED:
        raise CampaignError("job has unknown loss derivation")
    return loss_seed(job["loss_seed_base"], replicate)


def _empty_histogram():
    return {}


def _histogram_add(histogram, value, count=1):
    if (
        type(value) is not int or value < 0 or
        type(count) is not int or count < 0
    ):
        raise CampaignError("histogram input is invalid")
    key = str(value)
    histogram[key] = histogram.get(key, 0) + count


def _merge_histogram(destination, source):
    if not isinstance(source, dict):
        raise CampaignError("histogram is malformed")
    for key, count in source.items():
        if (
            type(key) is not str or not key.isdecimal() or
            str(int(key)) != key or type(count) is not int or count < 0
        ):
            raise CampaignError("histogram is malformed")
        destination[key] = destination.get(key, 0) + count


def _decode_compact_rows(columns, rows, expected_columns, expected_rows, label):
    if (
        columns != list(expected_columns) or
        not isinstance(rows, list) or len(rows) != expected_rows or
        any(not isinstance(row, list) or len(row) != len(columns)
            for row in rows)
    ):
        raise CampaignError(f"{label} compact rows are malformed")
    return [dict(zip(columns, row)) for row in rows]


def _new_compact_job_aggregate(job):
    parser_module = _require_v3_parser()
    stats_fields = tuple(parser_module.REPAIRTIMING_STATS_FIELDS)
    return {
        "schema": "wirehair.wh2.rv4a.compact-job-aggregate.v1",
        "job_id": job["job_id"],
        "phase": job["phase"],
        "lane": job["lane"],
        "arm": job["arm"],
        "arm_provisional_id": job["arm_provisional_id"],
        "K": job["K"],
        "schedule": job["schedule"],
        "cells": 0,
        "selector": {
            "observations": 0,
            "result_codes": {},
            "selected_attempt": {
                **{str(attempt): 0 for attempt in range(REPAIR_ATTEMPTS)},
                "none": 0,
            },
            "cap_exhausted": 0,
            "fatal_attempt_zero_mismatch": 0,
            "oom": 0,
            "uncommitted": 0,
            "execution_counts": {
                field: _empty_histogram()
                for field in SELECTOR_COUNT_FIELDS
            },
            "weak_constructions": {
                "raw_attempt0": 0,
                "repaired_final": 0,
            },
            "work": {
                field: {
                    "available": 0,
                    "unavailable": 0,
                    "histogram": _empty_histogram(),
                }
                for field in stats_fields
            },
        },
        "recovery": {
            role: {
                "outcome_classes": {},
                "success_overhead": [0] * (MAX_OVERHEAD + 1),
                "unrecovered": 0,
            }
            for role in ("raw", "repaired", "dispatch", "wh1")
        },
        "_role_hashers": {
            role: hashlib.sha256()
            for role in ("raw", "repaired", "dispatch", "wh1")
        },
        "candidate_final_weak": 0,
        "candidate_nonrepairable_regression": 0,
        "candidate_runtime_error": 0,
        "codec_oom": 0,
        "forced_control_failure": 0,
        "timing_accumulators": {
            panel: {
                "eligible_cells": 0,
                "censored_cells": 0,
                "mean": 0.0,
                "m2": 0.0,
            }
            for panel in parser_module.REPAIRTIMING_PANELS
        },
    }


def _accumulate_selector_work(aggregate, attempts, stats_fields):
    totals = {field: 0 for field in stats_fields}
    availability = {field: True for field in stats_fields}
    executed_calls = 0
    for expected_attempt, attempt in enumerate(attempts):
        if attempt.get("attempt") != expected_attempt:
            raise CampaignError("attempt rows are reordered")
        for kind in ("probe", "real"):
            executed = attempt.get(f"{kind}_executed")
            available = attempt.get(f"{kind}_stats_available")
            if executed not in (0, 1) or available not in (0, 1):
                raise CampaignError("attempt execution flags are invalid")
            if not executed:
                continue
            executed_calls += 1
            for field in stats_fields:
                value = attempt.get(f"{kind}_{field}")
                if available:
                    if type(value) is not int or value < 0:
                        raise CampaignError(
                            "available selector work is malformed")
                    totals[field] += value
                else:
                    availability[field] = False
    if executed_calls < 1:
        raise CampaignError("selector did not execute attempt-zero real work")
    for field in stats_fields:
        receipt = aggregate["selector"]["work"][field]
        if availability[field]:
            receipt["available"] += 1
            _histogram_add(receipt["histogram"], totals[field])
        else:
            receipt["unavailable"] += 1


def _timing_cell_values(real, parser_module):
    rows = _decode_compact_rows(
        real.get("timing_columns"),
        real.get("timing_rows"),
        parser_module.REPAIRTIMING_TIMING_FIELDS,
        len(parser_module.REPAIRTIMING_PANELS) * 8,
        "timing",
    )
    values = {}
    expected_order = "ABBABAAB"
    for panel_index, panel in enumerate(parser_module.REPAIRTIMING_PANELS):
        panel_rows = rows[panel_index * 8:(panel_index + 1) * 8]
        for expected_slot, row in enumerate(panel_rows):
            if (
                row.get("timing_panel") != panel or
                row.get("timing_panel_index") != panel_index or
                row.get("timing_slot") != expected_slot or
                row.get("timing_pair") != expected_slot // 2 or
                row.get("timing_label") != expected_order[expected_slot]
            ):
                raise CampaignError("timing rows are reordered")
        logs = []
        valid = True
        for pair in range(4):
            pair_rows = panel_rows[pair * 2:(pair + 1) * 2]
            by_label = {}
            for row in pair_rows:
                by_label[row["timing_label"]] = row
            if set(by_label) != {"A", "B"}:
                raise CampaignError("timing pair lacks A/B evidence")
            for row in by_label.values():
                if (
                    row.get("timing_eligible") != 1 or
                    row.get("timing_executed") != 1 or
                    row.get("timing_stable") != 1 or
                    row.get("timing_saturated") != 0 or
                    row.get("timing_cpu_migrated") != 0 or
                    row.get("timing_fault_contaminated") != 0 or
                    type(row.get("timing_elapsed_ns")) is not int or
                    row["timing_elapsed_ns"] < 1
                ):
                    valid = False
            if valid:
                logs.append(math.log(
                    by_label["A"]["timing_elapsed_ns"] /
                    by_label["B"]["timing_elapsed_ns"]))
        values[panel] = (
            math.fsum(logs) / 4.0 if valid and len(logs) == 4 else None)
    return values


def _accumulate_cell(
        compact, selector, real, *, accumulate_selector=True):
    parser_module = _require_v3_parser()
    selector_fields = {
        "schema", "arm_id", "K", "construction_root",
        "repair_policy_sha256", "selector_result", "selected_attempt",
        "attempts_executed", "calls_executed", "real_configuration_calls",
        "structural_probe_calls", "cap_exhausted",
        "fatal_attempt_zero_mismatch", "oom", "committed",
        "descriptor_sha256", "attempts", "selector_structural_sha256",
    }
    real_fields = {
        "schema", "trace_sha256", "message_sha256", "descriptor_hex",
        "roles", "controls", "attempt_columns", "attempt_rows",
        "timing_columns", "timing_rows",
    }
    if set(selector) != selector_fields or set(real) != real_fields:
        raise CampaignError("normalized v3 witness schema changed")
    if type(accumulate_selector) is not bool:
        raise CampaignError("selector aggregation flag is invalid")
    compact["cells"] += 1
    result = selector["selector_result"]
    if type(result) is not int or not 0 <= result <= 10:
        raise CampaignError("selector result code is invalid")
    if accumulate_selector:
        compact["selector"]["observations"] += 1
        _histogram_add(compact["selector"]["result_codes"], result)
    selected = selector["selected_attempt"]
    if selected == -1:
        if accumulate_selector:
            compact["selector"]["selected_attempt"]["none"] += 1
    elif type(selected) is int and 0 <= selected < REPAIR_ATTEMPTS:
        if accumulate_selector:
            compact["selector"]["selected_attempt"][str(selected)] += 1
    else:
        raise CampaignError("selected attempt is invalid")
    for source, destination in (
        ("cap_exhausted", "cap_exhausted"),
        ("fatal_attempt_zero_mismatch", "fatal_attempt_zero_mismatch"),
        ("oom", "oom"),
    ):
        value = selector[source]
        if value not in (0, 1):
            raise CampaignError("selector terminal flag is invalid")
        if accumulate_selector:
            compact["selector"][destination] += value
    committed = selector["committed"]
    if committed not in (0, 1):
        raise CampaignError("selector committed flag is invalid")
    if accumulate_selector:
        compact["selector"]["uncommitted"] += 1 - committed
    for field in SELECTOR_COUNT_FIELDS:
        value = selector[field]
        lower, upper = SELECTOR_COUNT_LIMITS[field]
        if type(value) is not int or not lower <= value <= upper:
            raise CampaignError(
                f"selector {field} accounting is invalid")
        if accumulate_selector:
            _histogram_add(
                compact["selector"]["execution_counts"][field], value)
    if (
        selector["calls_executed"] !=
            selector["real_configuration_calls"] +
            selector["structural_probe_calls"]
    ):
        raise CampaignError("selector call accounting is contradictory")

    attempts = _decode_compact_rows(
        real["attempt_columns"],
        real["attempt_rows"],
        parser_module.REPAIRTIMING_ATTEMPT_FIELDS,
        REPAIR_ATTEMPTS,
        "attempt",
    )
    if accumulate_selector:
        _accumulate_selector_work(
            compact, attempts, parser_module.REPAIRTIMING_STATS_FIELDS)

    roles = real["roles"]
    if not isinstance(roles, dict) or set(roles) != {
            "raw", "repaired", "dispatch", "wh1"}:
        raise CampaignError("recovery roles are malformed")
    for role, values in roles.items():
        if not isinstance(values, dict):
            raise CampaignError("recovery role is malformed")
        outcome = values.get("outcome_class")
        if type(outcome) is not str or not outcome:
            raise CampaignError("recovery outcome class is malformed")
        destination = compact["recovery"][role]
        destination["outcome_classes"][outcome] = (
            destination["outcome_classes"].get(outcome, 0) + 1)
        recovered = values.get("recovery_ok")
        overhead = values.get("overhead")
        if recovered == 1:
            if type(overhead) is not int or not 0 <= overhead <= MAX_OVERHEAD:
                raise CampaignError("recovery overhead is malformed")
            destination["success_overhead"][overhead] += 1
        elif recovered == 0 and overhead == -1:
            destination["unrecovered"] += 1
        else:
            raise CampaignError("recovery result is contradictory")
        compact["_role_hashers"][role].update(canonical_json_bytes({
            "cell": compact["cells"] - 1,
            "trace_sha256": real["trace_sha256"],
            "role": values,
        }))

    repaired = roles["repaired"]
    dispatch = roles["dispatch"]
    wh1 = roles["wh1"]
    repaired_outcome = repaired["outcome_class"]
    if accumulate_selector:
        compact["selector"]["weak_constructions"]["raw_attempt0"] += int(
            roles["raw"]["outcome_class"] == "weak")
        compact["selector"]["weak_constructions"]["repaired_final"] += int(
            repaired_outcome == "weak")
    if repaired_outcome == "weak":
        compact["candidate_final_weak"] += 1
    candidate_codes = (
        repaired.get("encode_result"),
        repaired.get("decode_construct_result"),
        repaired.get("feed_result"),
        repaired.get("recover_result"),
    )
    if (
        repaired.get("recovery_ok") != 1 and
        (
            dispatch.get("recovery_ok") == 1 or
            wh1.get("recovery_ok") == 1
        ) and
        repaired_outcome != "weak"
    ):
        compact["candidate_nonrepairable_regression"] += 1
    controls = real["controls"]
    if (
        not isinstance(controls, dict) or
        set(controls) != set(parser_module.REPAIRTIMING_CONTROL_FIELDS)
    ):
        raise CampaignError("forced/direct controls are malformed")
    execution_flags = (
        controls["repair_direct_executed"],
        controls["dispatch_direct_executed"],
    )
    if any(value not in (0, 1) for value in execution_flags):
        raise CampaignError("direct control execution flag is malformed")
    if (
        (committed and controls["forced_equal"] != 1) or
        (
            controls["repair_direct_executed"] == 1 and
            controls["repair_direct_witness_equal"] != 1
        ) or
        (
            controls["dispatch_direct_executed"] == 1 and
            controls["dispatch_direct_witness_equal"] != 1
        )
    ):
        compact["forced_control_failure"] += 1
    control_codes = tuple(
        controls[name] for name in (
            "forced_a_result", "forced_b_result",
            "repair_direct_result", "dispatch_direct_result",
        )
    )
    timing_rows = _decode_compact_rows(
        real["timing_columns"],
        real["timing_rows"],
        parser_module.REPAIRTIMING_TIMING_FIELDS,
        len(parser_module.REPAIRTIMING_PANELS) * 8,
        "timing",
    )
    timing_codes = tuple(
        row[name]
        for row in timing_rows
        for name in (
            "timing_construct_result",
            "timing_result",
            "timing_recover_result",
        )
    )
    role_codes = tuple(
        values.get(name)
        for values in roles.values()
        for name in (
            "encode_result", "decode_construct_result",
            "feed_result", "recover_result",
        )
    )
    attempt_codes = tuple(
        attempt[f"{kind}_result"]
        for attempt in attempts
        for kind in ("probe", "real")
        if attempt[f"{kind}_executed"] == 1
    )
    all_codes = (
        candidate_codes + control_codes + timing_codes +
        role_codes + attempt_codes
    )
    if 9 in all_codes:
        compact["codec_oom"] += 1
    if any(
            type(code) is int and
            code not in (-1, 0, 1, 3, 4, 7, 9)
            for code in all_codes):
        compact["candidate_runtime_error"] += 1
    for panel, value in _timing_cell_values(real, parser_module).items():
        accumulator = compact["timing_accumulators"][panel]
        if value is None:
            accumulator["censored_cells"] += 1
        else:
            _timing_accumulator_add(accumulator, value)


def _timing_accumulator_add(accumulator, value):
    """Add one finite value with Welford's cancellation-safe recurrence."""
    if type(value) is not float or not math.isfinite(value):
        raise CampaignError("timing accumulator value is invalid")
    count = accumulator["eligible_cells"]
    if type(count) is not int or count < 0:
        raise CampaignError("timing accumulator count is invalid")
    next_count = count + 1
    delta = value - accumulator["mean"]
    mean = accumulator["mean"] + delta / next_count
    accumulator["m2"] += delta * (value - mean)
    accumulator["mean"] = mean
    accumulator["eligible_cells"] = next_count


def _finish_timing_accumulator(accumulator):
    count = accumulator["eligible_cells"]
    if count < 1:
        return {
            **accumulator,
            "mean": None,
            "ci_low": None,
            "ci_high": None,
            "floor": None,
        }
    mean = accumulator["mean"]
    if count < 4:
        low = high = floor = None
    else:
        variance = max(
            0.0,
            accumulator["m2"] / (count - 1),
        )
        half = (
            peel_codec._student_t_critical_95(count - 1) *
            math.sqrt(variance / count)
        )
        low = mean - half
        high = mean + half
        floor = max(abs(low), abs(high))
    return {
        **accumulator,
        "mean": mean,
        "ci_low": low,
        "ci_high": high,
        "floor": floor,
    }


def _finish_compact_job_aggregate(compact):
    compact["recovery_role_stream_sha256"] = {
        role: digest.hexdigest()
        for role, digest in compact.pop("_role_hashers").items()
    }
    compact["timing"] = {
        panel: _finish_timing_accumulator(accumulator)
        for panel, accumulator in compact.pop(
            "timing_accumulators").items()
    }


def _cell_aggregate(parser_module, evidence, job):
    """Reduce a receipt to queue-bounded hashes and sufficient statistics."""
    cells = getattr(evidence, "cells", None)
    if not isinstance(cells, tuple) or len(cells) != job["replicates"]:
        raise CampaignError("repairtiming evidence cell count is invalid")
    selectors = {}
    real_digest = hashlib.sha256()
    compact = _new_compact_job_aggregate(job)
    for replicate, cell in enumerate(cells):
        if (
            type(getattr(cell, "replicate", None)) is not int or
            cell.replicate != replicate or
            type(getattr(cell, "construction_root", None)) is not int or
            cell.construction_root != _expected_root(job, replicate) or
            type(getattr(cell, "loss_seed", None)) is not int or
            cell.loss_seed != _expected_loss(job, replicate)
        ):
            raise CampaignError("repairtiming cell seed/order changed")
        expected_key = (
            job["arm_provisional_id"], job["K"], cell.construction_root)
        if getattr(cell, "selector_key", None) != expected_key:
            raise CampaignError("repairtiming selector key changed")
        selector = parser_module.selector_projection(cell)
        real = parser_module.real_projection(cell)
        if not isinstance(selector, dict) or not isinstance(real, dict):
            raise CampaignError("repairtiming projection is malformed")
        if (
            selector.get("schema") !=
                "wirehair.wh2.repair-v1.selector-structural.v1" or
            selector.get("arm_id") !=
                int(job["arm_provisional_id"], 16) or
            selector.get("K") != job["K"] or
            selector.get("construction_root") != cell.construction_root or
            selector.get("repair_policy_sha256") !=
                REPAIR_POLICY_SHA256 or
            real.get("schema") !=
                "wirehair.wh2.repairtiming.real-witness.v3"
        ):
            raise CampaignError(
                "repairtiming projection changed its frozen coordinates")
        prior = selectors.get(expected_key)
        is_new_selector = prior is None
        if prior is None:
            selectors[expected_key] = selector
        elif not peel_codec._same_typed_json(prior, selector):
            raise CampaignError(
                "duplicate selector key has payload-dependent structure")
        real_digest.update(canonical_json_bytes({
            "replicate": replicate,
            "construction_root": cell.construction_root,
            "loss_seed": cell.loss_seed,
            "real": real,
        }))
        _accumulate_cell(
            compact,
            selector,
            real,
            accumulate_selector=is_new_selector,
        )
    selector_rows = [
        {
            "selector_key": list(key),
            "selector_projection": selectors[key],
        }
        for key in sorted(selectors)
    ]
    compact["selector_set_sha256"] = canonical_sha256(selector_rows)
    compact["unique_selectors"] = len(selector_rows)
    if compact["selector"]["observations"] != compact["unique_selectors"]:
        raise CampaignError("compact selector deduplication lost a key")
    compact["real_stream_sha256"] = real_digest.hexdigest()
    _finish_compact_job_aggregate(compact)
    return compact


def _receipt_context_projection(receipt, assignment):
    context = receipt.get("context")
    if not isinstance(context, dict):
        raise CampaignError("receipt context is missing")
    bound = context.get("bound")
    thermal = context.get("thermal")
    if (
        not isinstance(bound, dict) or
        bound.get("cpu_affinity") != [assignment["cpu"]] or
        not isinstance(thermal, dict) or
        type(thermal.get("valid_rows")) is not int or
        thermal["valid_rows"] < 2 or
        thermal.get("edac_ce_max") != 0 or
        thermal.get("edac_ue_max") != 0
    ):
        raise CampaignError("receipt thermal/EDAC coverage is invalid")
    return {
        "context_sha256": canonical_sha256(context),
        "cpu": assignment["cpu"],
        "thermal_device": bound.get("thermal_device"),
        "thermal_inode": bound.get("thermal_inode"),
        "valid_rows": thermal["valid_rows"],
        "cpu_tctl_max_c": thermal.get("cpu_tctl_max_c"),
        "dimm_max_c": thermal.get("dimm_max_c"),
        "edac_ce_max": 0,
        "edac_ue_max": 0,
        "first_monotonic_s": thermal.get("first_monotonic_s"),
        "last_monotonic_s": thermal.get("last_monotonic_s"),
    }


def _replay_result(directory, manifest, manifest_sha256, assignment, job):
    parser_module = _require_v3_parser()
    limits = _result_evidence_byte_limits(job)
    path = directory / assignment["output"]
    envelope, evidence = read_canonical_gzip_json(
        path,
        compressed_limit=limits["compressed"],
        uncompressed_limit=limits["uncompressed"],
    )
    expected_fields = {
        "schema", "manifest_sha256", "job", "assignment",
        "wall_started_unix_ns", "wall_finished_unix_ns",
        "runtime_bindings_before", "runtime_bindings_after", "receipt",
    }
    if (
        not isinstance(envelope, dict) or set(envelope) != expected_fields or
        envelope.get("schema") != RESULT_SCHEMA or
        envelope.get("manifest_sha256") != manifest_sha256 or
        not peel_codec._same_typed_json(envelope.get("job"), job) or
        not peel_codec._same_typed_json(
            envelope.get("assignment"), assignment) or
        type(envelope.get("wall_started_unix_ns")) is not int or
        type(envelope.get("wall_finished_unix_ns")) is not int or
        not 0 < envelope["wall_started_unix_ns"] <= \
            envelope["wall_finished_unix_ns"] or
        not peel_codec._same_typed_json(
            envelope.get("runtime_bindings_before"),
            manifest["runtime_bindings"]) or
        not peel_codec._same_typed_json(
            envelope.get("runtime_bindings_after"),
            manifest["runtime_bindings"])
    ):
        raise CampaignError(f"result envelope is invalid: {path}")
    receipt = envelope["receipt"]
    replayed = parser_module.replay_repairtiming_receipt(
        receipt, expected_request=expected_request(job))
    if (
        getattr(replayed, "protocol", None) != REQUIRED_REPAIRTIMING_PROTOCOL or
        not peel_codec._same_typed_json(
            _measurement_as_dict(replayed), receipt)
    ):
        raise CampaignError("strict result replay changed receipt")
    thermal = _receipt_context_projection(receipt, assignment)
    if (
        thermal["thermal_device"] !=
            manifest["runtime_bindings"]["thermal"]["device"] or
        thermal["thermal_inode"] !=
            manifest["runtime_bindings"]["thermal"]["inode"]
    ):
        raise CampaignError("result receipt changed thermal identity")
    aggregate = _cell_aggregate(parser_module, replayed, job)
    return {
        "order": assignment["order"],
        "job_id": job["job_id"],
        "path": assignment["output"],
        **evidence,
        "wall_started_unix_ns": envelope["wall_started_unix_ns"],
        "wall_finished_unix_ns": envelope["wall_finished_unix_ns"],
        "thermal": thermal,
        "aggregate": aggregate,
    }


_REPLAY_PROCESS_STATE = None


def _initialize_replay_process(
        directory, manifest, manifest_sha256):
    global _REPLAY_PROCESS_STATE
    for item in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        signal.signal(item, signal.SIG_DFL)
    _REPLAY_PROCESS_STATE = (
        Path(directory), manifest, manifest_sha256)


def _replay_process_task(task):
    if _REPLAY_PROCESS_STATE is None:
        raise CampaignError("replay worker was not initialized")
    index, assignment, job = task
    directory, manifest, manifest_sha256 = _REPLAY_PROCESS_STATE
    result = _replay_result(
        directory, manifest, manifest_sha256, assignment, job)
    # The process queue carries only one compact aggregate and hashes, never
    # raw native rows, normalized cells, or attempt lists.
    return index, result


def _replay_process_main(
        directory, manifest, manifest_sha256, tasks, results):
    _initialize_replay_process(directory, manifest, manifest_sha256)
    try:
        while True:
            task = tasks.get()
            if task is None:
                return
            index, result = _replay_process_task(task)
            results.put({"kind": "result", "index": index, "value": result})
    except BaseException as error:
        try:
            results.put({
                "kind": "error",
                "pid": os.getpid(),
                "type": type(error).__name__,
                "message": str(error),
            })
        except BaseException:
            pass
        raise


def _stop_replay_processes(processes):
    deadline = time.monotonic() + REPLAY_TERMINATE_GRACE_SECONDS
    for process in processes:
        if process.pid is not None and process.is_alive():
            process.terminate()
    for process in processes:
        if process.pid is None:
            continue
        remaining = max(0.0, deadline - time.monotonic())
        process.join(remaining)
    deadline = time.monotonic() + REPLAY_KILL_GRACE_SECONDS
    for process in processes:
        if process.pid is not None and process.is_alive():
            process.kill()
    for process in processes:
        if process.pid is None:
            continue
        remaining = max(0.0, deadline - time.monotonic())
        process.join(remaining)
    survivors = [
        process.pid for process in processes
        if process.pid is not None and process.is_alive()
    ]
    if survivors:
        raise CampaignError(
            f"parallel replay workers survived SIGKILL: {survivors!r}")


def _ordered_replayed_results(
        directory, manifest, manifest_sha256, replay_workers):
    assignments = manifest["assignments"]
    jobs = manifest["pre_cpu_jobs"]
    count = len(jobs)
    if (
        type(replay_workers) is not int or
        not 1 <= replay_workers <= MAX_WORKERS
    ):
        raise CampaignError("replay worker count is invalid")
    worker_count = min(replay_workers, count)
    context = multiprocessing.get_context("fork")
    tasks = context.Queue(maxsize=REPLAY_QUEUE_MAX_ITEMS)
    results = context.Queue(maxsize=REPLAY_QUEUE_MAX_ITEMS)
    processes = [
        context.Process(
            target=_replay_process_main,
            args=(
                str(directory), manifest, manifest_sha256, tasks, results),
            daemon=False,
        )
        for unused_index in range(worker_count)
    ]
    pending = {}
    submitted = received = yielded = 0
    try:
        for process in processes:
            process.start()
        while yielded < count:
            _parent_signal_safe_point()
            while (
                submitted < count and
                submitted - yielded < REPLAY_QUEUE_MAX_ITEMS
            ):
                try:
                    tasks.put_nowait((
                        submitted,
                        assignments[submitted],
                        jobs[submitted],
                    ))
                except queue.Full:
                    break
                submitted += 1
            try:
                message = results.get(timeout=REPLAY_QUEUE_POLL_SECONDS)
            except queue.Empty:
                failed = [
                    (process.pid, process.exitcode)
                    for process in processes
                    if process.exitcode is not None
                ]
                if failed:
                    raise CampaignError(
                        "parallel replay worker exited before shutdown: "
                        f"{failed!r}")
                continue
            if (
                not isinstance(message, dict) or
                message.get("kind") not in ("result", "error")
            ):
                raise CampaignError("parallel replay queue was corrupted")
            if message["kind"] == "error":
                raise CampaignError(
                    "parallel replay worker error "
                    f"{message.get('type')}: {message.get('message')}")
            index = message.get("index")
            if (
                type(index) is not int or not 0 <= index < count or
                index in pending or index < yielded
            ):
                raise CampaignError("parallel replay result order is invalid")
            pending[index] = message["value"]
            received += 1
            if (
                len(pending) > REPLAY_QUEUE_MAX_ITEMS or
                submitted - yielded > REPLAY_QUEUE_MAX_ITEMS
            ):
                raise CampaignError(
                    "parallel replay reorder window exceeded its bound")
            while yielded in pending:
                yield pending.pop(yielded)
                yielded += 1
        for unused_index in processes:
            while True:
                _parent_signal_safe_point()
                try:
                    tasks.put(
                        None, timeout=REPLAY_QUEUE_POLL_SECONDS)
                    break
                except queue.Full:
                    failed = [
                        (process.pid, process.exitcode)
                        for process in processes
                        if process.exitcode not in (None, 0)
                    ]
                    if failed:
                        raise CampaignError(
                            "parallel replay worker failed during "
                            f"shutdown: {failed!r}")
        deadline = time.monotonic() + REPLAY_SHUTDOWN_GRACE_SECONDS
        while any(process.is_alive() for process in processes):
            _parent_signal_safe_point()
            for process in processes:
                process.join(0)
            failed = [
                (process.pid, process.exitcode)
                for process in processes
                if process.exitcode not in (None, 0)
            ]
            if failed:
                raise CampaignError(
                    f"parallel replay worker failed: {failed!r}")
            if time.monotonic() >= deadline:
                raise CampaignError(
                    "parallel replay workers did not exit after sentinels")
            time.sleep(REPLAY_QUEUE_POLL_SECONDS)
        failed = [
            (process.pid, process.exitcode)
            for process in processes
            if process.exitcode != 0
        ]
        if failed:
            raise CampaignError(
                f"parallel replay worker failed: {failed!r}")
    except BaseException:
        _stop_replay_processes(processes)
        raise
    finally:
        for queue_object in (tasks, results):
            try:
                queue_object.cancel_join_thread()
                queue_object.close()
            except BaseException:
                pass


class AtomicGzipJsonLines:
    """Streaming, deterministic JSONL publication with no replacement."""

    def __init__(
            self, path, *, row_limit, compressed_limit,
            uncompressed_limit):
        self.path = Path(path)
        self.temporary = self.path.with_name(self.path.name + ".tmp")
        if (
            type(row_limit) is not int or row_limit < 0 or
            type(compressed_limit) is not int or compressed_limit < 1 or
            type(uncompressed_limit) is not int or uncompressed_limit < 1
        ):
            raise CampaignError("JSONL publication limits are invalid")
        self.row_limit = row_limit
        self.compressed_limit = compressed_limit
        self.uncompressed_limit = uncompressed_limit
        self.raw = None
        self.archive = None
        self.compressed_hash = hashlib.sha256()
        self.uncompressed_hash = hashlib.sha256()
        self.compressed_bytes = 0
        self.uncompressed_bytes = 0
        self.rows = 0

    def __enter__(self):
        flags = os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
        descriptor = os.open(str(self.temporary), flags, 0o666)
        self.raw = os.fdopen(descriptor, "w+b")
        self.archive = gzip.GzipFile(
            filename="", mode="wb", fileobj=self.raw,
            compresslevel=6, mtime=0)
        return self

    def write(self, value):
        payload = canonical_json_bytes(value)
        if self.rows >= self.row_limit:
            raise CampaignError("JSONL publication exceeds row cap")
        if self.uncompressed_bytes + len(payload) > self.uncompressed_limit:
            raise CampaignError(
                "JSONL publication exceeds uncompressed cap")
        self.archive.write(payload)
        self.uncompressed_hash.update(payload)
        self.uncompressed_bytes += len(payload)
        self.rows += 1

    def __exit__(self, error_type, error, traceback):
        try:
            try:
                if self.archive is not None:
                    self.archive.close()
            finally:
                if self.raw is not None and not self.raw.closed:
                    try:
                        self.raw.flush()
                        os.fsync(self.raw.fileno())
                        if error_type is None:
                            stable = (
                                "st_dev", "st_ino", "st_mode", "st_size",
                                "st_mtime_ns", "st_ctime_ns",
                            )
                            before = os.fstat(self.raw.fileno())
                            if (
                                not stat.S_ISREG(before.st_mode) or
                                before.st_size > self.compressed_limit
                            ):
                                raise CampaignError(
                                    "JSONL publication exceeds compressed "
                                    "cap")
                            self.raw.seek(0)
                            while True:
                                chunk = self.raw.read(1024 * 1024)
                                if not chunk:
                                    break
                                self.compressed_hash.update(chunk)
                                self.compressed_bytes += len(chunk)
                            after = os.fstat(self.raw.fileno())
                            current = os.stat(
                                self.temporary, follow_symlinks=False)
                            if (
                                any(getattr(before, name) !=
                                    getattr(after, name)
                                    for name in stable) or
                                any(getattr(after, name) !=
                                    getattr(current, name)
                                    for name in stable) or
                                self.compressed_bytes != after.st_size
                            ):
                                raise CampaignError(
                                    "JSONL publication changed before "
                                    "commit")
                    finally:
                        self.raw.close()
            if error_type is not None:
                try:
                    os.unlink(self.temporary)
                except FileNotFoundError:
                    pass
                return False
            _commit_no_replace(self.temporary, self.path)
        except BaseException:
            try:
                os.unlink(self.temporary)
            except FileNotFoundError:
                pass
            raise
        return False

    def evidence(self, relative_path):
        return {
            "path": relative_path,
            "rows": self.rows,
            "compressed_sha256": self.compressed_hash.hexdigest(),
            "uncompressed_sha256": self.uncompressed_hash.hexdigest(),
            "compressed_bytes": self.compressed_bytes,
            "uncompressed_bytes": self.uncompressed_bytes,
        }


def _new_phase_bucket(name):
    parser_module = _require_v3_parser()
    return {
        "name": name,
        "jobs": 0,
        "cells": 0,
        "recovery": {
            role: {
                "outcome_classes": {},
                "success_overhead": [0] * (MAX_OVERHEAD + 1),
                "unrecovered": 0,
            }
            for role in ("raw", "repaired", "dispatch", "wh1")
        },
        "selector": {
            "observations": 0,
            "result_codes": {},
            "selected_attempt": {
                **{str(attempt): 0 for attempt in range(REPAIR_ATTEMPTS)},
                "none": 0,
            },
            "cap_exhausted": 0,
            "fatal_attempt_zero_mismatch": 0,
            "oom": 0,
            "uncommitted": 0,
            "execution_counts": {
                field: {}
                for field in SELECTOR_COUNT_FIELDS
            },
            "weak_constructions": {
                "raw_attempt0": 0,
                "repaired_final": 0,
            },
            "work": {
                field: {
                    "available": 0,
                    "unavailable": 0,
                    "histogram": {},
                }
                for field in parser_module.REPAIRTIMING_STATS_FIELDS
            },
        },
        "candidate_final_weak": 0,
        "candidate_nonrepairable_regression": 0,
        "candidate_runtime_error": 0,
        "codec_oom": 0,
        "forced_control_failure": 0,
        "timing_jobs": [],
        "thermal": {
            "valid_rows": 0,
            "cpu_tctl_max_c": 0.0,
            "dimm_max_c": 0.0,
            "edac_ce_max": 0,
            "edac_ue_max": 0,
        },
    }


def _merge_counts(destination, source):
    if not isinstance(source, dict):
        raise CampaignError("count map is malformed")
    for key, value in source.items():
        if type(key) is not str or type(value) is not int or value < 0:
            raise CampaignError("count map is malformed")
        destination[key] = destination.get(key, 0) + value


def _merge_compact_into_bucket(
        bucket, compact, thermal, *, include_selector=True):
    if type(include_selector) is not bool:
        raise CampaignError("selector merge flag is invalid")
    bucket["jobs"] += 1
    bucket["cells"] += compact["cells"]
    for role in bucket["recovery"]:
        target = bucket["recovery"][role]
        source = compact["recovery"][role]
        _merge_counts(
            target["outcome_classes"], source["outcome_classes"])
        if (
            not isinstance(source["success_overhead"], list) or
            len(source["success_overhead"]) != MAX_OVERHEAD + 1 or
            any(type(value) is not int or value < 0
                for value in source["success_overhead"]) or
            type(source["unrecovered"]) is not int or
            source["unrecovered"] < 0
        ):
            raise CampaignError("recovery histogram is malformed")
        target["success_overhead"] = [
            left + right for left, right in zip(
                target["success_overhead"],
                source["success_overhead"],
            )
        ]
        target["unrecovered"] += source["unrecovered"]
    if include_selector:
        selector_target = bucket["selector"]
        selector_source = compact["selector"]
        observations = selector_source.get("observations")
        if (
            type(observations) is not int or observations < 0 or
            observations != compact["unique_selectors"]
        ):
            raise CampaignError(
                "compact selector observation count is malformed")
        selector_target["observations"] += observations
        _merge_histogram(
            selector_target["result_codes"],
            selector_source["result_codes"],
        )
        for key, value in selector_source["selected_attempt"].items():
            if key not in selector_target["selected_attempt"] or \
                    type(value) is not int or value < 0:
                raise CampaignError(
                    "selected attempt histogram is malformed")
            selector_target["selected_attempt"][key] += value
        for field in (
                "cap_exhausted", "fatal_attempt_zero_mismatch", "oom",
                "uncommitted"):
            selector_target[field] += selector_source[field]
        if (
            set(selector_source.get("execution_counts", {})) !=
                set(SELECTOR_COUNT_FIELDS) or
            set(selector_source.get("weak_constructions", {})) !=
                {"raw_attempt0", "repaired_final"}
        ):
            raise CampaignError(
                "selector count or weak-construction evidence is malformed")
        for field in SELECTOR_COUNT_FIELDS:
            _merge_histogram(
                selector_target["execution_counts"][field],
                selector_source["execution_counts"][field],
            )
        for field in ("raw_attempt0", "repaired_final"):
            value = selector_source["weak_constructions"][field]
            if type(value) is not int or not 0 <= value <= observations:
                raise CampaignError(
                    "weak-construction evidence is malformed")
            selector_target["weak_constructions"][field] += value
        for field, source in selector_source["work"].items():
            target = selector_target["work"][field]
            target["available"] += source["available"]
            target["unavailable"] += source["unavailable"]
            _merge_histogram(target["histogram"], source["histogram"])
    for field in (
        "candidate_final_weak", "candidate_nonrepairable_regression",
        "candidate_runtime_error", "codec_oom", "forced_control_failure",
    ):
        bucket[field] += compact[field]
    bucket["timing_jobs"].append({
        "job_id": compact["job_id"],
        "K": compact["K"],
        "schedule": compact["schedule"],
        "lane": compact["lane"],
        "cells": compact["cells"],
        "panels": compact["timing"],
    })
    for field in ("valid_rows",):
        if type(thermal.get(field)) is not int or thermal[field] < 0:
            raise CampaignError("thermal aggregate is malformed")
        bucket["thermal"][field] += thermal[field]
    for field in ("cpu_tctl_max_c", "dimm_max_c"):
        value = thermal.get(field)
        if (
            type(value) not in (int, float) or
            isinstance(value, bool) or not math.isfinite(value)
        ):
            raise CampaignError("thermal aggregate is malformed")
        bucket["thermal"][field] = max(
            bucket["thermal"][field], float(value))
    for field in ("edac_ce_max", "edac_ue_max"):
        if thermal.get(field) != 0:
            raise CampaignError("EDAC error found in authenticated evidence")


def _histogram_statistics(histogram, available, unavailable):
    if (
        type(available) is not int or available < 0 or
        type(unavailable) is not int or unavailable < 0
    ):
        raise CampaignError("work availability is malformed")
    ordered = []
    total = 0
    weighted = 0
    for key in sorted(histogram, key=int):
        count = histogram[key]
        value = int(key)
        if type(count) is not int or count < 0:
            raise CampaignError("work histogram is malformed")
        ordered.append((value, count))
        total += count
        weighted += value * count
    if total != available:
        raise CampaignError("work histogram lost available observations")
    if available == 0:
        return {
            "available": 0,
            "unavailable": unavailable,
            "mean": None,
            "p50": None,
            "p95": None,
            "p99": None,
            "max": None,
        }

    def nearest_rank(probability):
        rank = max(1, math.ceil(probability * available))
        cumulative = 0
        for value, count in ordered:
            cumulative += count
            if cumulative >= rank:
                return value
        raise CampaignError("work percentile histogram is incomplete")

    return {
        "available": available,
        "unavailable": unavailable,
        "mean": weighted / available,
        "p50": nearest_rank(0.50),
        "p95": nearest_rank(0.95),
        "p99": nearest_rank(0.99),
        "max": next(
            value for value, count in reversed(ordered) if count > 0),
    }


def _selected_attempt_statistics(histogram):
    selected_count = sum(
        histogram[str(attempt)] for attempt in range(REPAIR_ATTEMPTS))
    expanded = {
        str(attempt): histogram[str(attempt)]
        for attempt in range(REPAIR_ATTEMPTS)
    }
    statistics = _histogram_statistics(
        expanded, selected_count, histogram["none"])
    return {
        "histogram": dict(histogram),
        "selected_statistics": statistics,
    }


def _validate_selector_aggregate(selector):
    observations = selector.get("observations")
    if (
        type(observations) is not int or observations < 1 or
        sum(selector["result_codes"].values()) != observations or
        sum(selector["selected_attempt"].values()) != observations or
        any(
            type(selector[field]) is not int or
            not 0 <= selector[field] <= observations
            for field in (
                "cap_exhausted", "fatal_attempt_zero_mismatch",
                "oom", "uncommitted",
            )
        ) or
        any(
            receipt["available"] + receipt["unavailable"] != observations
            for receipt in selector["work"].values()
        ) or
        set(selector.get("execution_counts", {})) !=
            set(SELECTOR_COUNT_FIELDS) or
        any(
            not isinstance(histogram, dict) or
            any(
                type(key) is not str or not key.isdecimal() or
                len(key) > 2 or str(int(key)) != key or
                not lower <= int(key) <= upper or
                type(count) is not int or count < 0
                for key, count in histogram.items()
            ) or
            sum(histogram.values()) != observations
            for field, histogram in selector["execution_counts"].items()
            for lower, upper in (SELECTOR_COUNT_LIMITS[field],)
        ) or
        set(selector.get("weak_constructions", {})) !=
            {"raw_attempt0", "repaired_final"} or
        any(
            type(value) is not int or not 0 <= value <= observations
            for value in selector["weak_constructions"].values()
        )
    ):
        raise CampaignError(
            "selector evidence does not cover every unique selector")


def _recovery_metrics(receipt, cells):
    if (
        not isinstance(receipt, dict) or
        set(receipt) != {
            "outcome_classes", "success_overhead", "unrecovered"
        }
    ):
        raise CampaignError("recovery evidence schema is malformed")
    outcome_classes = receipt["outcome_classes"]
    success_overhead = receipt["success_overhead"]
    unrecovered = receipt["unrecovered"]
    if (
        type(cells) is not int or cells < 1 or
        not isinstance(outcome_classes, dict) or
        any(
            type(name) is not str or not name or
            type(count) is not int or count < 0
            for name, count in outcome_classes.items()
        ) or
        sum(outcome_classes.values()) != cells or
        not isinstance(success_overhead, list) or
        len(success_overhead) != MAX_OVERHEAD + 1 or
        any(type(count) is not int or count < 0
            for count in success_overhead) or
        type(unrecovered) is not int or unrecovered < 0 or
        sum(success_overhead) + unrecovered != cells
    ):
        raise CampaignError("recovery evidence does not cover every cell")
    cumulative = 0
    tail_auc = 0
    for overhead, count in enumerate(success_overhead):
        if overhead > MAX_OVERHEAD:
            raise CampaignError("recovery overhead exceeded cap")
        cumulative += count
        tail_auc += cells - cumulative
    return {
        "cells": cells,
        "outcome_classes": dict(outcome_classes),
        "success_overhead": list(success_overhead),
        "final_unrecovered_at_h64": unrecovered,
        "overhead_tail_auc_h0_to_h64": tail_auc,
    }


def _student_t_interval_from_values(values):
    if (
        not isinstance(values, list) or len(values) < 4 or
        any(type(value) is not float or not math.isfinite(value)
            for value in values)
    ):
        raise CampaignError("timing aggregate lacks finite job means")
    count = len(values)
    mean = math.fsum(values) / count
    variance = math.fsum(
        (value - mean) ** 2 for value in values) / (count - 1)
    half = (
        peel_codec._student_t_critical_95(count - 1) *
        math.sqrt(variance / count)
    )
    return mean, mean - half, mean + half


TIMING_CONFIRMATIONS = {
    "full_encoder_candidate_wh1": (
        "encoder_selector_wh1",
        "encoder_selector_aa",
        "encoder_wh1_aa",
    ),
    "decoder_feed_candidate_wh1": (
        "decoder_feed_wh1",
        "decoder_feed_aa",
        "decoder_feed_wh1_aa",
    ),
    "full_decoder_candidate_wh1": (
        "decoder_full_wh1",
        "decoder_full_aa",
        "decoder_full_wh1_aa",
    ),
    "selected_direct_candidate_dispatch": (
        "direct_dispatch",
        "direct_aa",
        "direct_dispatch_aa",
    ),
}
TIMING_REPORTS = {
    **TIMING_CONFIRMATIONS,
    "full_encoder_selector_forced": (
        "encoder_selector_forced",
        "encoder_selector_aa",
        "encoder_forced_aa",
    ),
}


def _timing_confirmation(jobs, metric):
    cross_name, left_aa_name, right_aa_name = TIMING_REPORTS[metric]
    required_lower, unused_upper = \
        peel_codec._paired_practical_log_bounds(REQUIRED_MARGIN)
    required_floor = -required_lower
    means = []
    authenticated_highs = []
    effective_floors = []
    incomplete = []
    for job in jobs:
        panels = job["panels"]
        try:
            cross = panels[cross_name]
            left_aa = panels[left_aa_name]
            right_aa = panels[right_aa_name]
        except KeyError:
            raise CampaignError("timing panel roster changed")
        minimum = job["cells"] // 2
        available = all(
            panel["eligible_cells"] >= minimum and
            type(panel["mean"]) is float and math.isfinite(panel["mean"]) and
            type(panel["ci_low"]) is float and
            math.isfinite(panel["ci_low"]) and
            type(panel["ci_high"]) is float and
            math.isfinite(panel["ci_high"]) and
            type(panel["floor"]) is float and math.isfinite(panel["floor"])
            for panel in (cross, left_aa, right_aa)
        )
        if not available:
            incomplete.append(job["job_id"])
            continue
        means.append(cross["mean"])
        authenticated_highs.append(cross["ci_high"])
        effective_floors.append(max(
            required_floor, left_aa["floor"], right_aa["floor"]))
    if incomplete:
        return {
            "complete": False,
            "jobs": len(jobs),
            "complete_jobs": len(means),
            "incomplete_job_ids": incomplete,
            "mean": None,
            "ci_low": None,
            "ci_high": None,
            "mean_authenticated_job_ci_high": None,
            "aggregate_upper": None,
            "mean_effective_floor": None,
            "pass": False,
        }
    mean, low, high = _student_t_interval_from_values(means)
    mean_high = math.fsum(authenticated_highs) / len(authenticated_highs)
    aggregate_upper = max(high, mean_high)
    mean_floor = math.fsum(effective_floors) / len(effective_floors)
    return {
        "complete": True,
        "jobs": len(jobs),
        "complete_jobs": len(means),
        "incomplete_job_ids": [],
        "mean": mean,
        "ci_low": low,
        "ci_high": high,
        "mean_authenticated_job_ci_high": mean_high,
        "aggregate_upper": aggregate_upper,
        "mean_effective_floor": mean_floor,
        "pass": aggregate_upper < -mean_floor,
    }


def _finalize_phase_bucket(bucket):
    _validate_selector_aggregate(bucket["selector"])
    recovery = {
        role: _recovery_metrics(receipt, bucket["cells"])
        for role, receipt in bucket["recovery"].items()
    }
    selector = {
        key: value for key, value in bucket["selector"].items()
        if key not in ("work", "execution_counts")
    }
    selector["selected_attempt"] = _selected_attempt_statistics(
        bucket["selector"]["selected_attempt"])
    selector["work"] = {
        field: _histogram_statistics(
            receipt["histogram"],
            receipt["available"],
            receipt["unavailable"],
        )
        for field, receipt in bucket["selector"]["work"].items()
    }
    selector["execution_counts"] = {
        field: {
            "histogram": dict(histogram),
            "statistics": _histogram_statistics(
                histogram,
                selector["observations"],
                0,
            ),
        }
        for field, histogram in
        bucket["selector"]["execution_counts"].items()
    }
    timing = {
        metric: _timing_confirmation(bucket["timing_jobs"], metric)
        for metric in TIMING_REPORTS
    }
    return {
        "name": bucket["name"],
        "jobs": bucket["jobs"],
        "cells": bucket["cells"],
        "recovery": recovery,
        "selector": selector,
        "candidate_final_weak": bucket["candidate_final_weak"],
        "candidate_nonrepairable_regression":
            bucket["candidate_nonrepairable_regression"],
        "candidate_runtime_error": bucket["candidate_runtime_error"],
        "codec_oom": bucket["codec_oom"],
        "forced_control_failure": bucket["forced_control_failure"],
        "timing": timing,
        "thermal": dict(bucket["thermal"]),
    }


def _recovery_no_worse(candidate, control):
    return (
        candidate["final_unrecovered_at_h64"] <=
            control["final_unrecovered_at_h64"] and
        candidate["overhead_tail_auc_h0_to_h64"] <=
            control["overhead_tail_auc_h0_to_h64"]
    )


def _arm_eligibility(summary):
    recovery = summary["recovery"]
    selector = summary["selector"]
    zero_counts = {
        "selector_cap_exhaustion": selector["cap_exhausted"],
        "selector_fatal_attempt_zero_mismatch":
            selector["fatal_attempt_zero_mismatch"],
        "selector_uncommitted": selector["uncommitted"],
        "candidate_final_weak": summary["candidate_final_weak"],
        "candidate_nonrepairable_regression":
            summary["candidate_nonrepairable_regression"],
        "candidate_runtime_error": summary["candidate_runtime_error"],
        "selector_or_codec_OOM": selector["oom"] + summary["codec_oom"],
        "forced_control_failure": summary["forced_control_failure"],
    }
    recovery_pass = (
        _recovery_no_worse(recovery["repaired"], recovery["dispatch"]) and
        _recovery_no_worse(recovery["repaired"], recovery["wh1"])
    )
    timing_pass = {
        metric: summary["timing"][metric]["pass"]
        for metric in TIMING_CONFIRMATIONS
    }
    return {
        "zero_counts": zero_counts,
        "zero_gate_pass": all(value == 0 for value in zero_counts.values()),
        "recovery_pass": recovery_pass,
        "timing_pass": timing_pass,
        "eligible": (
            all(value == 0 for value in zero_counts.values()) and
            recovery_pass and all(timing_pass.values())
        ),
    }


def _derive_sealed_confirmation(survivor, overall, lane_summaries):
    """Apply frozen gates to the union and both lanes without adaptation."""
    if survivor not in ARM_BY_NAME or (
        not isinstance(lane_summaries, dict) or
        set(lane_summaries) != {"random", "production"}
    ):
        raise CampaignError("sealed confirmation inputs are malformed")
    eligibility = _arm_eligibility(overall)
    lane_eligibility = {
        lane: _arm_eligibility(lane_summaries[lane])
        for lane in ("random", "production")
    }
    confirmed = (
        eligibility["eligible"] and
        all(item["eligible"] for item in lane_eligibility.values())
    )
    return {
        "schema": "wirehair.wh2.rv4a.sealed-confirmation.v1",
        "survivor": survivor,
        "eligibility": eligibility,
        "lane_eligibility": lane_eligibility,
        "status": "confirmed" if confirmed else "failed",
        "lane_breakdowns": lane_summaries,
        "public_promotion": "not-evaluated-by-rv4a2",
        "fallback": "forbidden",
    }


def _derive_training_decision(arm_summaries):
    candidates = {}
    eligible = []
    for arm_name, summary in arm_summaries.items():
        eligibility = _arm_eligibility(summary)
        candidates[arm_name] = {
            "eligibility": eligibility,
            "fresh_attempt0_raw": summary["recovery"]["raw"],
            "fresh_weak_constructions":
                dict(summary["selector"]["weak_constructions"]),
            "historical_za5v_raw_weak_constructions":
                ARM_BY_NAME[arm_name][
                    "historical_za5v_raw_weak_constructions"],
            "repaired": summary["recovery"]["repaired"],
            "dispatch": summary["recovery"]["dispatch"],
            "wh1": summary["recovery"]["wh1"],
            "timing": summary["timing"],
        }
        if eligibility["eligible"]:
            eligible.append(arm_name)
    eligible.sort(key=lambda arm_name: (
        arm_summaries[arm_name]["timing"][
            "decoder_feed_candidate_wh1"]["mean"],
        arm_summaries[arm_name]["timing"][
            "full_decoder_candidate_wh1"]["mean"],
        arm_summaries[arm_name]["timing"][
            "selected_direct_candidate_dispatch"]["mean"],
        arm_summaries[arm_name]["timing"][
            "full_encoder_candidate_wh1"]["mean"],
        arm_name,
    ))
    selected = eligible[0] if eligible else None
    return {
        "schema": TRAINING_DECISION_SCHEMA,
        "policy_sha256": canonical_sha256(selection_policy()),
        "candidate_order": eligible,
        "candidates": candidates,
        "status": "winner" if selected is not None else "no-survivor",
        "selected_survivor": selected,
        "sealed_authorization":
            "one-joint-sealed-run" if selected is not None else "forbidden",
        "public_promotion": "not-evaluated-by-rv4a2",
    }


def _validate_cross_job_invariants(compact, selector_groups, controls):
    selector_key = (
        compact["arm"], compact["lane"], compact["K"])
    current = (
        compact["selector_set_sha256"], compact["unique_selectors"])
    prior = selector_groups.get(selector_key)
    if prior is None:
        selector_groups[selector_key] = {
            "value": current,
            "schedules": {compact["schedule"]},
        }
    else:
        if prior["value"] != current:
            raise CampaignError(
                "selector projection changed across loss schedules")
        if compact["schedule"] in prior["schedules"]:
            raise CampaignError("duplicate selector schedule job")
        prior["schedules"].add(compact["schedule"])
    if compact["phase"] == "training":
        control_key = (compact["lane"], compact["K"], compact["schedule"])
        value = (
            compact["recovery_role_stream_sha256"]["dispatch"],
            compact["recovery_role_stream_sha256"]["wh1"],
        )
        prior_control = controls.get(control_key)
        if prior_control is None:
            controls[control_key] = value
        elif prior_control != value:
            raise CampaignError(
                "candidate-independent recovery controls drifted across arms")


def _validate_selector_groups(selector_groups, phase):
    unique = 0
    for key, receipt in selector_groups.items():
        if receipt["schedules"] != set(SCHEDULES):
            raise CampaignError("selector group lacks all six schedules")
        unique += receipt["value"][1]
    expected = (
        152064 if phase == "training" else 25443)
    if unique != expected:
        raise CampaignError(
            f"unique selector cardinality changed: {unique} != {expected}")
    return unique


def build_summaries(
        directory, manifest, output_directory, *, replay_workers=MAX_WORKERS):
    """Replay every result and atomically build compact authenticated evidence."""
    _require_source_forced_outcome_runner()
    directory = Path(directory).resolve()
    output_directory = Path(output_directory).resolve()
    manifest_sha256 = canonical_sha256(manifest)
    jobs = _validate_manifest(manifest)
    _require_source_forced_outcome_runner(
        manifest["runtime_bindings"]["sources"]["runner"]["sha256"])
    job_files = _verify_job_files(
        directory, manifest, manifest_sha256)
    result_files = []
    selector_groups = {}
    controls = {}
    arm_buckets = {
        arm_name: _new_phase_bucket(f"arm:{arm_name}")
        for arm_name in (
            ARM_BY_NAME if manifest["phase"] == "training"
            else (manifest["survivor"],)
        )
    }
    lane_buckets = (
        {
            lane: _new_phase_bucket(f"lane:{lane}")
            for lane in ("random", "production")
        }
        if manifest["phase"] == "sealed" else {}
    )
    overall = _new_phase_bucket("overall")
    ledger_path = output_directory / "cell-ledger.jsonl.gz"
    with AtomicGzipJsonLines(
            ledger_path,
            row_limit=len(jobs),
            compressed_limit=SUMMARY_COMPRESSED_BYTE_LIMIT,
            uncompressed_limit=SUMMARY_UNCOMPRESSED_BYTE_LIMIT,
    ) as ledger:
        for result in _ordered_replayed_results(
                directory, manifest, manifest_sha256, replay_workers):
            _parent_signal_safe_point()
            compact = result["aggregate"]
            _validate_cross_job_invariants(
                compact, selector_groups, controls)
            # Selector metrics use exactly one deterministic real-work
            # witness per frozen (arm,K,root): the first replicate carrying
            # that key in the first frozen schedule.  Recovery and timing
            # remain weighted per physical cell.
            include_selector = compact["schedule"] == SCHEDULES[0]
            _merge_compact_into_bucket(
                arm_buckets[compact["arm"]],
                compact,
                result["thermal"],
                include_selector=include_selector,
            )
            _merge_compact_into_bucket(
                overall,
                compact,
                result["thermal"],
                include_selector=include_selector,
            )
            if lane_buckets:
                _merge_compact_into_bucket(
                    lane_buckets[compact["lane"]],
                    compact,
                    result["thermal"],
                    include_selector=include_selector,
                )
            result_files.append({
                key: value for key, value in result.items()
                if key != "aggregate"
            })
            ledger.write({
                "schema": "wirehair.wh2.rv4a.job-ledger.v1",
                "order": result["order"],
                "job_id": result["job_id"],
                "result_compressed_sha256":
                    result["compressed_sha256"],
                "result_uncompressed_sha256":
                    result["uncompressed_sha256"],
                "selector_set_sha256":
                    compact["selector_set_sha256"],
                "unique_selectors": compact["unique_selectors"],
                "real_stream_sha256": compact["real_stream_sha256"],
                "recovery_role_stream_sha256":
                    compact["recovery_role_stream_sha256"],
                "compact_aggregate_sha256":
                    canonical_sha256(compact),
                "thermal_context_sha256":
                    result["thermal"]["context_sha256"],
            })
    unique_selectors = _validate_selector_groups(
        selector_groups, manifest["phase"])
    if overall["selector"]["observations"] != unique_selectors:
        raise CampaignError(
            "selector metrics are not weighted by unique selector")
    finalized_arms = {
        name: _finalize_phase_bucket(bucket)
        for name, bucket in arm_buckets.items()
    }
    finalized_overall = _finalize_phase_bucket(overall)
    finalized_lanes = {
        name: _finalize_phase_bucket(bucket)
        for name, bucket in lane_buckets.items()
    }
    if manifest["phase"] == "training":
        decision = _derive_training_decision(finalized_arms)
        sealed_confirmation = None
    else:
        decision = None
        sealed_confirmation = _derive_sealed_confirmation(
            manifest["survivor"],
            finalized_overall,
            finalized_lanes,
        )
    summary = {
        "schema": SUMMARY_SCHEMA,
        "phase": manifest["phase"],
        "survivor": manifest["survivor"],
        "manifest_sha256": manifest_sha256,
        "frozen_roster_sha256": manifest["frozen_roster_sha256"],
        "pre_cpu_job_list_sha256":
            manifest["pre_cpu_job_list_sha256"],
        "cell_set_sha256": manifest["cell_set_sha256"],
        "policy_sha256": manifest["policy_sha256"],
        "preflight_pin_sha256": manifest["preflight_pin_sha256"],
        "jobs": len(jobs),
        "cells": sum(job["replicates"] for job in jobs),
        "unique_selectors": unique_selectors,
        "attempt_rows":
            sum(job["replicates"] for job in jobs) * REPAIR_ATTEMPTS,
        "native_rows":
            sum(job["replicates"] for job in jobs) * 121,
        "arms": finalized_arms,
        "overall": finalized_overall,
        "training_decision": decision,
        "sealed_confirmation": sealed_confirmation,
    }
    summary_path = output_directory / "summary.json.gz"
    summary_evidence = atomic_gzip_json(summary_path, summary, limits={
        "compressed": SUMMARY_COMPRESSED_BYTE_LIMIT,
        "uncompressed": SUMMARY_UNCOMPRESSED_BYTE_LIMIT,
    })
    return {
        "job_files": job_files,
        "result_files": result_files,
        "cell_ledger": ledger.evidence("cell-ledger.jsonl.gz"),
        "summary": {
            "path": "summary.json.gz",
            **summary_evidence,
        },
        "training_decision_sha256":
            None if decision is None else canonical_sha256(decision),
        "sealed_confirmation_sha256":
            None if sealed_confirmation is None
            else canonical_sha256(sealed_confirmation),
    }


def _atomic_bytes(path, payload):
    path = Path(path)
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


def _write_completion_checksum(directory, digest):
    if not _is_sha256(digest):
        raise CampaignError("completion digest is invalid")
    _atomic_bytes(
        directory / "completion.sha256",
        f"{digest}  completion.json\n".encode("ascii"),
    )


def _validate_training_contract_shape(contract, *, expected_survivor):
    fields = {
        "path", "manifest_sha256", "completion_sha256",
        "training_decision_sha256", "selected_survivor",
        "sealed_job_list_sha256", "sealed_cell_set_sha256",
    }
    validation = validate_frozen_contract()
    if (
        not isinstance(contract, dict) or set(contract) != fields or
        not isinstance(contract.get("path"), str) or
        not Path(contract["path"]).is_absolute() or
        any(not _is_sha256(contract.get(field)) for field in (
            "manifest_sha256", "completion_sha256",
            "training_decision_sha256", "sealed_job_list_sha256",
            "sealed_cell_set_sha256",
        )) or
        contract.get("selected_survivor") != expected_survivor or
        expected_survivor not in ARM_BY_NAME or
        contract["sealed_job_list_sha256"] !=
            validation["sealed_job_list_sha256"][expected_survivor] or
        contract["sealed_cell_set_sha256"] !=
            validation["sealed_cell_set_sha256"][expected_survivor]
    ):
        raise CampaignError("sealed training contract is invalid")
    return contract


def _load_summary(directory, completion):
    path = Path(directory) / "summary.json.gz"
    summary, evidence = read_canonical_gzip_json(
        path,
        compressed_limit=SUMMARY_COMPRESSED_BYTE_LIMIT,
        uncompressed_limit=SUMMARY_UNCOMPRESSED_BYTE_LIMIT,
    )
    expected = completion.get("summary")
    if (
        not isinstance(expected, dict) or
        not peel_codec._same_typed_json(expected, {
            "path": "summary.json.gz",
            **evidence,
        }) or
        not isinstance(summary, dict) or
        summary.get("schema") != SUMMARY_SCHEMA
    ):
        raise CampaignError("stored summary does not match completion")
    return summary


def _load_training_contract(directory):
    directory = Path(directory).resolve()
    verified = verify_campaign(directory)
    manifest, manifest_sha256 = read_canonical_json(
        directory / "manifest.json")
    completion, completion_sha256 = read_canonical_json(
        directory / "completion.json")
    summary = _load_summary(directory, completion)
    decision = summary.get("training_decision")
    survivor = (
        decision.get("selected_survivor")
        if isinstance(decision, dict) else None
    )
    if (
        manifest.get("phase") != "training" or
        completion.get("phase") != "training" or
        verified.get("phase") != "training" or
        not isinstance(decision, dict) or
        decision.get("schema") != TRAINING_DECISION_SCHEMA or
        decision.get("status") != "winner" or
        decision.get("sealed_authorization") != "one-joint-sealed-run" or
        survivor not in ARM_BY_NAME or
        completion.get("training_decision_sha256") !=
            canonical_sha256(decision)
    ):
        raise CampaignError(
            "training evidence does not authorize a sealed campaign")
    validation = validate_frozen_contract()
    contract = {
        "path": str(directory),
        "manifest_sha256": manifest_sha256,
        "completion_sha256": completion_sha256,
        "training_decision_sha256": canonical_sha256(decision),
        "selected_survivor": survivor,
        "sealed_job_list_sha256":
            validation["sealed_job_list_sha256"][survivor],
        "sealed_cell_set_sha256":
            validation["sealed_cell_set_sha256"][survivor],
    }
    return _validate_training_contract_shape(
        contract, expected_survivor=survivor)


def _validate_runtime_monitor(monitor, manifest):
    fields = {
        "scheduler", "checks", "checks_while_running",
        "first_check_monotonic_ns", "last_check_monotonic_ns",
        "runtime_bindings_sha256", "startup_slots",
        "pending_idle_checks", "drain_policy",
        "drained_slots_while_global_unstarted_checks",
        "completed_groups_confirmed_gone",
        "occupancy_histogram", "cpu_queues",
    }
    if (
        not isinstance(monitor, dict) or
        set(monitor) != fields or
        monitor.get("scheduler") != "per-cpu-rolling-v1" or
        type(monitor.get("checks")) is not int or
        monitor["checks"] < 2 or
        type(monitor.get("checks_while_running")) is not int or
        not 1 <= monitor["checks_while_running"] < monitor["checks"] or
        type(monitor.get("first_check_monotonic_ns")) is not int or
        type(monitor.get("last_check_monotonic_ns")) is not int or
        not 0 < monitor["first_check_monotonic_ns"] <=
            monitor["last_check_monotonic_ns"] or
        monitor.get("runtime_bindings_sha256") !=
            canonical_sha256(manifest["runtime_bindings"]) or
        monitor.get("startup_slots") != manifest["workers"] or
        monitor.get("pending_idle_checks") != 0 or
        monitor.get("drain_policy") !=
            "idle-only-after-that-cpu-fixed-queue-exhausted-v1" or
        type(monitor.get(
            "drained_slots_while_global_unstarted_checks")) is not int or
        monitor["drained_slots_while_global_unstarted_checks"] < 0 or
        monitor["drained_slots_while_global_unstarted_checks"] >
            monitor["checks_while_running"] or
        monitor.get("completed_groups_confirmed_gone") !=
            len(manifest["assignments"]) or
        not isinstance(monitor.get("occupancy_histogram"), dict) or
        not isinstance(monitor.get("cpu_queues"), list) or
        len(monitor["cpu_queues"]) != manifest["workers"]
    ):
        raise CampaignError("runtime monitor is malformed")

    occupancy_total = occupancy_running = 0
    for key, count in monitor["occupancy_histogram"].items():
        if (
            type(key) is not str or not key.isdecimal() or
            str(int(key)) != key or
            not 0 <= int(key) <= manifest["workers"] or
            type(count) is not int or count < 1
        ):
            raise CampaignError(
                "runtime monitor occupancy histogram is malformed")
        occupancy_total += count
        if int(key) > 0:
            occupancy_running += count
    if (
        occupancy_total != monitor["checks"] or
        occupancy_running != monitor["checks_while_running"] or
        monitor["occupancy_histogram"].get("0") != 1 or
        monitor["occupancy_histogram"].get(
            str(manifest["workers"]), 0) < 1
    ):
        raise CampaignError("runtime monitor occupancy totals are invalid")

    by_cpu = {
        cpu: [
            assignment["order"]
            for assignment in manifest["assignments"]
            if assignment["cpu"] == cpu
        ]
        for cpu in manifest["worker_cpus"]
    }
    seen = []
    earliest_start = None
    latest_finish = None
    for expected_cpu, queue_receipt in zip(
            manifest["worker_cpus"], monitor["cpu_queues"]):
        if (
            not isinstance(queue_receipt, dict) or
            set(queue_receipt) != {"cpu", "orders", "launches"} or
            queue_receipt.get("cpu") != expected_cpu or
            queue_receipt.get("orders") != by_cpu[expected_cpu] or
            not isinstance(queue_receipt.get("launches"), list) or
            len(queue_receipt["launches"]) != len(by_cpu[expected_cpu])
        ):
            raise CampaignError("runtime monitor CPU queue is malformed")
        prior_finish = None
        for expected_order, launch in zip(
                by_cpu[expected_cpu], queue_receipt["launches"]):
            if (
                not isinstance(launch, dict) or set(launch) != {
                    "order", "started_monotonic_ns",
                    "finished_monotonic_ns", "returncode",
                } or
                launch.get("order") != expected_order or
                type(launch.get("started_monotonic_ns")) is not int or
                type(launch.get("finished_monotonic_ns")) is not int or
                not 0 < launch["started_monotonic_ns"] <=
                    launch["finished_monotonic_ns"] or
                launch.get("returncode") != 0 or
                (
                    prior_finish is not None and
                    launch["started_monotonic_ns"] < prior_finish
                )
            ):
                raise CampaignError(
                    "runtime monitor worker launch is malformed")
            prior_finish = launch["finished_monotonic_ns"]
            earliest_start = (
                launch["started_monotonic_ns"]
                if earliest_start is None else
                min(earliest_start, launch["started_monotonic_ns"])
            )
            latest_finish = (
                launch["finished_monotonic_ns"]
                if latest_finish is None else
                max(latest_finish, launch["finished_monotonic_ns"])
            )
            seen.append(expected_order)
    if (
        sorted(seen) != list(range(len(manifest["assignments"]))) or
        earliest_start is None or latest_finish is None or
        monitor["first_check_monotonic_ns"] < earliest_start or
        monitor["last_check_monotonic_ns"] < latest_finish
    ):
        raise CampaignError("runtime monitor lost worker launches")


def _run_campaign_impl(
        phase, directory, bench_path, thermal_source, *,
        pin_record, pin_sha256, training_directory=None,
        workers=MAX_WORKERS, replay_workers=MAX_WORKERS):
    _require_source_forced_outcome_runner()
    _require_normal_priority()
    validate_frozen_contract()
    if (
        type(workers) is not int or not 1 <= workers <= MAX_WORKERS or
        type(replay_workers) is not int or
        not 1 <= replay_workers <= MAX_WORKERS
    ):
        raise CampaignError("worker count is invalid")
    if phase == "training":
        if training_directory is not None:
            raise CampaignError("training forbids parent evidence")
        training_contract = None
        survivor = None
    elif phase == "sealed":
        if training_directory is None:
            raise CampaignError(
                "sealed requires authenticated training evidence")
        training_contract = _load_training_contract(training_directory)
        survivor = training_contract["selected_survivor"]
    else:
        raise CampaignError(f"unknown phase {phase!r}")
    cpus = sorted(os.sched_getaffinity(0))
    if len(cpus) < workers:
        raise CampaignError(
            f"requested {workers} workers but only {len(cpus)} CPUs allowed")
    cpus = cpus[:workers]
    jobs = build_pre_cpu_jobs(phase, survivor)
    bindings = runtime_bindings(bench_path, thermal_source)
    authenticated_pin_sha256 = _validate_preflight_pin(
        pin_record, pin_sha256, bindings)
    directory = _create_fresh_directory(Path(directory).resolve())
    for name in ("jobs", "results", "logs"):
        (directory / name).mkdir()
    manifest = _build_manifest(
        phase,
        survivor,
        jobs,
        cpus,
        bindings,
        training_contract,
        authenticated_pin_sha256,
    )
    manifest_sha256 = atomic_json(directory / "manifest.json", manifest)
    job_files = _write_job_files(
        directory, manifest, manifest_sha256)
    _parent_signal_safe_point()
    check_runtime_bindings(bindings, full_hash=True)
    monitor = _run_rolling_workers(directory, manifest)
    check_runtime_bindings(bindings, full_hash=True)
    _parent_signal_safe_point()
    logs = _verify_empty_logs(directory, manifest)
    summaries = build_summaries(
        directory,
        manifest,
        directory,
        replay_workers=replay_workers,
    )
    if not peel_codec._same_typed_json(summaries["job_files"], job_files):
        raise CampaignError("job files changed during campaign")
    final_bindings = runtime_bindings(
        bindings["benchmark"]["path"], bindings["thermal"]["path"])
    if not peel_codec._same_typed_json(final_bindings, bindings):
        raise CampaignError("runtime binding changed before completion")
    completion = {
        "schema": COMPLETION_SCHEMA,
        "phase": phase,
        "survivor": survivor,
        "manifest_sha256": manifest_sha256,
        "frozen_roster_sha256": manifest["frozen_roster_sha256"],
        "pre_cpu_job_list_sha256":
            manifest["pre_cpu_job_list_sha256"],
        "cell_set_sha256": manifest["cell_set_sha256"],
        "policy_sha256": manifest["policy_sha256"],
        "preflight_pin_sha256": manifest["preflight_pin_sha256"],
        "jobs": len(jobs),
        "runtime_binding_monitor": monitor,
        "runtime_bindings_final": final_bindings,
        "job_files": summaries["job_files"],
        "log_files": logs,
        "result_files": summaries["result_files"],
        "cell_ledger": summaries["cell_ledger"],
        "summary": summaries["summary"],
        "training_decision_sha256":
            summaries["training_decision_sha256"],
        "sealed_confirmation_sha256":
            summaries["sealed_confirmation_sha256"],
        "completed_unix_ns": time.time_ns(),
    }
    completion_sha256 = atomic_json(
        directory / "completion.json", completion)
    _write_completion_checksum(directory, completion_sha256)
    _parent_signal_safe_point()
    return completion


def run_campaign(
        phase, directory, bench_path, thermal_source, *,
        pin_record, pin_sha256, training_directory=None,
        workers=MAX_WORKERS, replay_workers=MAX_WORKERS):
    with _parent_termination_handlers():
        completion = _run_campaign_impl(
            phase,
            directory,
            bench_path,
            thermal_source,
            pin_record=pin_record,
            pin_sha256=pin_sha256,
            training_directory=training_directory,
            workers=workers,
            replay_workers=replay_workers,
        )
        _parent_signal_safe_point()
        return completion


def _stream_file_sha256(path, *, byte_limit):
    payload = _read_stable_bytes(path, byte_limit=byte_limit)
    return sha256_bytes(payload), len(payload)


def _verify_stored_artifacts(directory, completion):
    summary = _load_summary(directory, completion)
    ledger = completion.get("cell_ledger")
    if (
        not isinstance(ledger, dict) or set(ledger) != {
            "path", "rows", "compressed_sha256", "uncompressed_sha256",
            "compressed_bytes", "uncompressed_bytes",
        } or
        ledger.get("path") != "cell-ledger.jsonl.gz" or
        type(ledger.get("rows")) is not int or
        ledger["rows"] != completion["jobs"] or
        type(ledger.get("compressed_bytes")) is not int or
        not 0 <= ledger["compressed_bytes"] <=
            SUMMARY_COMPRESSED_BYTE_LIMIT or
        type(ledger.get("uncompressed_bytes")) is not int or
        not 0 <= ledger["uncompressed_bytes"] <=
            SUMMARY_UNCOMPRESSED_BYTE_LIMIT or
        not _is_sha256(ledger.get("compressed_sha256")) or
        not _is_sha256(ledger.get("uncompressed_sha256"))
    ):
        raise CampaignError("stored ledger evidence is malformed")
    rows, replayed = read_canonical_gzip_json_lines(
        directory / ledger["path"],
        compressed_limit=SUMMARY_COMPRESSED_BYTE_LIMIT,
        uncompressed_limit=SUMMARY_UNCOMPRESSED_BYTE_LIMIT,
        row_limit=completion["jobs"],
    )
    if (
        len(rows) != completion["jobs"] or
        not peel_codec._same_typed_json(replayed, {
            "compressed_sha256": ledger["compressed_sha256"],
            "uncompressed_sha256": ledger["uncompressed_sha256"],
            "compressed_bytes": ledger["compressed_bytes"],
            "uncompressed_bytes": ledger["uncompressed_bytes"],
            "rows": ledger["rows"],
        })
    ):
        raise CampaignError("stored ledger does not match completion")
    return summary


def _verify_campaign_impl(directory, *, replay_workers):
    _require_source_forced_outcome_runner()
    directory = Path(directory).resolve()
    manifest, manifest_sha256 = read_canonical_json(
        directory / "manifest.json")
    _validate_manifest(manifest)
    _require_source_forced_outcome_runner(
        manifest["runtime_bindings"]["sources"]["runner"]["sha256"])
    if manifest["phase"] == "sealed":
        reopened = _load_training_contract(
            manifest["training_contract"]["path"])
        if not peel_codec._same_typed_json(
                reopened, manifest["training_contract"]):
            raise CampaignError("sealed parent training evidence changed")
    completion, completion_sha256 = read_canonical_json(
        directory / "completion.json")
    checksum_payload = _read_stable_bytes(
        directory / "completion.sha256", byte_limit=256)
    try:
        checksum = checksum_payload.decode("ascii")
    except UnicodeDecodeError as error:
        raise CampaignError(
            f"completion checksum is not ASCII: {error}")
    if checksum != f"{completion_sha256}  completion.json\n":
        raise CampaignError("completion checksum does not match")
    fields = {
        "schema", "phase", "survivor", "manifest_sha256",
        "frozen_roster_sha256", "pre_cpu_job_list_sha256",
        "cell_set_sha256", "policy_sha256", "preflight_pin_sha256", "jobs",
        "runtime_binding_monitor", "runtime_bindings_final",
        "job_files", "log_files", "result_files", "cell_ledger", "summary",
        "training_decision_sha256", "sealed_confirmation_sha256",
        "completed_unix_ns",
    }
    if (
        not isinstance(completion, dict) or set(completion) != fields or
        completion.get("schema") != COMPLETION_SCHEMA or
        completion.get("phase") != manifest["phase"] or
        completion.get("survivor") != manifest["survivor"] or
        completion.get("manifest_sha256") != manifest_sha256 or
        completion.get("frozen_roster_sha256") !=
            manifest["frozen_roster_sha256"] or
        completion.get("pre_cpu_job_list_sha256") !=
            manifest["pre_cpu_job_list_sha256"] or
        completion.get("cell_set_sha256") != manifest["cell_set_sha256"] or
        completion.get("policy_sha256") != manifest["policy_sha256"] or
        completion.get("preflight_pin_sha256") !=
            manifest["preflight_pin_sha256"] or
        completion.get("jobs") != len(manifest["pre_cpu_jobs"]) or
        not peel_codec._same_typed_json(
            completion.get("runtime_bindings_final"),
            manifest["runtime_bindings"]) or
        type(completion.get("completed_unix_ns")) is not int or
        completion["completed_unix_ns"] < manifest["created_unix_ns"]
    ):
        raise CampaignError("completion does not match manifest")
    _validate_runtime_monitor(
        completion["runtime_binding_monitor"], manifest)
    stored_summary = _verify_stored_artifacts(directory, completion)
    check_runtime_bindings(manifest["runtime_bindings"], full_hash=True)
    with tempfile.TemporaryDirectory(
            prefix="wh2-rv4a-verify-") as temporary:
        reproduced = build_summaries(
            directory,
            manifest,
            Path(temporary),
            replay_workers=replay_workers,
        )
    for field in (
        "job_files", "result_files", "cell_ledger", "summary",
        "training_decision_sha256", "sealed_confirmation_sha256",
    ):
        if not peel_codec._same_typed_json(
                reproduced[field], completion[field]):
            raise CampaignError(
                f"strict replay did not reproduce {field}")
    if not peel_codec._same_typed_json(
            _verify_empty_logs(directory, manifest),
            completion["log_files"]):
        raise CampaignError("worker logs changed after completion")
    decision = stored_summary.get("training_decision")
    confirmation = stored_summary.get("sealed_confirmation")
    if (
        completion["training_decision_sha256"] != (
            None if decision is None else canonical_sha256(decision)) or
        completion["sealed_confirmation_sha256"] != (
            None if confirmation is None
            else canonical_sha256(confirmation))
    ):
        raise CampaignError("stored decision hash is invalid")
    return {
        "phase": manifest["phase"],
        "survivor": manifest["survivor"],
        "manifest_sha256": manifest_sha256,
        "completion_sha256": completion_sha256,
        "jobs": completion["jobs"],
        "cells": sum(
            job["replicates"] for job in manifest["pre_cpu_jobs"]),
        "training_decision_sha256":
            completion["training_decision_sha256"],
        "sealed_confirmation_sha256":
            completion["sealed_confirmation_sha256"],
    }


def verify_campaign(directory, *, replay_workers=MAX_WORKERS):
    _require_source_forced_outcome_runner()
    if _PARENT_SIGNAL_STATE is not None:
        return _verify_campaign_impl(
            directory, replay_workers=replay_workers)
    with _parent_termination_handlers():
        result = _verify_campaign_impl(
            directory, replay_workers=replay_workers)
        _parent_signal_safe_point()
        return result


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser(
        "plan", help="print a result-free frozen plan")
    preflight = subparsers.add_parser(
        "preflight-pin",
        help="print complete source/binary pins before any outcomes",
    )
    preflight.add_argument("--bench", required=True, type=Path)
    preflight.add_argument(
        "--thermal-source", required=True, type=Path)
    training = subparsers.add_parser(
        "run-training", help="execute the frozen two-arm training phase")
    training.add_argument("--directory", required=True, type=Path)
    training.add_argument("--bench", required=True, type=Path)
    training.add_argument("--thermal-source", required=True, type=Path)
    training.add_argument("--pin-record", required=True, type=Path)
    training.add_argument("--pin-sha256", required=True)
    training.add_argument("--workers", type=int, default=MAX_WORKERS)
    training.add_argument(
        "--replay-workers", type=int, default=MAX_WORKERS)

    sealed = subparsers.add_parser(
        "run-sealed",
        help="execute one decision-derived joint sealed boundary",
    )
    sealed.add_argument("--directory", required=True, type=Path)
    sealed.add_argument("--bench", required=True, type=Path)
    sealed.add_argument("--thermal-source", required=True, type=Path)
    sealed.add_argument("--pin-record", required=True, type=Path)
    sealed.add_argument("--pin-sha256", required=True)
    sealed.add_argument(
        "--training-directory", required=True, type=Path)
    sealed.add_argument("--workers", type=int, default=MAX_WORKERS)
    sealed.add_argument(
        "--replay-workers", type=int, default=MAX_WORKERS)

    verify = subparsers.add_parser(
        "verify", help="strictly replay one completed boundary")
    verify.add_argument("--directory", required=True, type=Path)
    verify.add_argument(
        "--replay-workers", type=int, default=MAX_WORKERS)

    worker = subparsers.add_parser("_worker")
    worker.add_argument("--job-file", required=True, type=Path)
    worker.add_argument("--output", required=True, type=Path)
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    if args.command == "plan":
        print(canonical_json_bytes(
            make_result_free_plan()
        ).decode("ascii"), end="")
    elif args.command == "preflight-pin":
        print(canonical_json_bytes(
            make_preflight_pin(
                args.bench, args.thermal_source)
        ).decode("ascii"), end="")
    elif args.command == "run-training":
        run_campaign(
            "training",
            args.directory,
            args.bench,
            args.thermal_source,
            pin_record=args.pin_record,
            pin_sha256=args.pin_sha256,
            workers=args.workers,
            replay_workers=args.replay_workers,
        )
    elif args.command == "run-sealed":
        run_campaign(
            "sealed",
            args.directory,
            args.bench,
            args.thermal_source,
            pin_record=args.pin_record,
            pin_sha256=args.pin_sha256,
            training_directory=args.training_directory,
            workers=args.workers,
            replay_workers=args.replay_workers,
        )
    elif args.command == "verify":
        print(canonical_json_bytes(
            verify_campaign(
                args.directory,
                replay_workers=args.replay_workers,
            )
        ).decode("ascii"), end="")
    else:
        run_worker(args.job_file, args.output)


def _source_forced_cli(argv):
    """Run CLI logic from one exact self-source snapshot with a hash stamp."""
    path = Path(__file__).resolve()
    payload = _bootstrap_source_bytes(
        path, byte_limit=PYTHON_SOURCE_BYTE_LIMIT)
    digest = sha256_bytes(payload)
    module_name = f"_wirehair_rv4a_cli_{digest[:16]}"
    module = types.ModuleType(module_name)
    module.__file__ = str(path)
    module.__package__ = ""
    module.__loader__ = None
    module.__spec__ = None
    module._BOOTSTRAP_RUNNER_SOURCE_SHA256 = digest
    sys.modules[module_name] = module
    try:
        code = compile(payload, str(path), "exec", dont_inherit=True)
        exec(code, module.__dict__)
        if module._bootstrap_source_bytes(
                path,
                byte_limit=module.PYTHON_SOURCE_BYTE_LIMIT,
        ) != payload:
            raise RuntimeError(
                "campaign runner changed during source-forced launch")
        module.main(argv)
    except Exception as error:
        known = tuple(
            item for item in (
                getattr(module, "CampaignError", None),
                getattr(
                    getattr(module, "peel_codec", None),
                    "MeasurementError",
                    None,
                ),
                ValueError,
            )
            if isinstance(item, type)
        )
        if isinstance(error, known):
            print(f"error: {error}", file=sys.stderr)
            return 2
        raise
    finally:
        sys.modules.pop(module_name, None)
    return 0


if __name__ == "__main__":
    raise SystemExit(_source_forced_cli(sys.argv[1:]))
