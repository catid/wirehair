#!/usr/bin/env python3
"""Freeze and run post-selection H12 q0 preferred-attempt research.

This harness never changes a production profile.  It evaluates a broad seed
route: for each K, try a preferred joint precode/packet ``seed_attempt`` at the
actual block width; an invalid attempt falls back to the canonical ascending
first-valid attempt and a preferred attempt equal to the canonical attempt is
an exact no-op alias.  Only systematic-valid, non-no-op direct observations
can win development ranking.

The complete protocol, source cohort, timing-derived 128-bin assignment,
external loss roots, round panels, ranking tuple, H1 pruning gates, H2
whole-table gates, and unsaturated timing contract are sealed before R1.  A
future control root is cached physically but masked from every earlier
selection.  H1 may only prune a routed K back to canonical control; it may not
substitute a runner-up.  H2 and timing may only accept or reject the whole
frozen table.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import fcntl
import hashlib
import heapq
import importlib
import io
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import stat as stat_module
import subprocess
import sys
import tempfile
import threading
import time
from typing import Any, Iterable, Iterator, Sequence
import unicodedata

import wh2_rank_floor_two_anchor_allk as common
from wh2_rank_floor_two_anchor_screen import (
    THERMAL_FIELDS,
    thermal_finish,
    thermal_start,
)


SCHEMA = "wirehair.wh2.h12_preferred_attempt.contract.v1"
SOURCE_COMMIT = "df13c42966233ac64ae6b65e4aed318485313ad8"
SOURCE_FAILURES_SHA256 = (
    "3f0a28c33eea4a707ac38136ef9ff5125d4fb5502655022734ab680b549dca39"
)
SOURCE_PAIRED_SHA256 = (
    "ab4a1f1e31b59a302f9803c99bed29bb1bfbb1715a2a38d470196c5356d41587"
)
SOURCE_PHASE_SEAL_SHA256 = (
    "8c30c922849e2e13ff0c3e696415071ea9985f42c372642ce62cc50fd4092823"
)
SOURCE_ANALYSIS_SEAL_SHA256 = (
    "98bd961ad2badac4d4c3fad344644546a015c8054164aa6b21ac9e8ab83960f2"
)
SOURCE_CONTRACT_SHA256 = (
    "5f9e48809a4656e90e236edabbe7bd3c20bbfc4ee85ef222c6fae1fc4c8bc9ca"
)
SOURCE_ROOTS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
)
COHORT_SHA256 = (
    "1233920b35cd2b594ec48ef7dcd016a0c5132f28b79a2c914b7d843989fea6c5"
)
WEIGHTS_SHA256 = (
    "ab99a2a39de5b5c22ac6e4ce30885e16ed9c29d40449961668e856558a6e9c4b"
)
ASSIGNMENT_SHA256 = (
    "1450eefbd5b1cdc2542d6de5b08dd1beed48f18ece35fd6927b052c8f44343f7"
)
BINS_SHA256 = (
    "e55dd00bd63b81a919975b965b7bae8ebd00edae4a80260f20d7c2ac74106a83"
)
ALLK_ASSIGNMENT_SHA256 = (
    "d483f00f9d5564f3a3e7638c0aaa00ff2b8500494e56fe2299014fb21023388d"
)
ALLK_BINS_SHA256 = (
    "6ee99c364654bc020959b4ae5f20edeae703a54615fb3de57530b8e9b2a25430"
)

K_MIN = 2
K_MAX = 64000
CUTOFF = 4096
ACTIVE_K_COUNT = 641
BIN_COUNT = 128
WIDTHS = (64, 256, 1280, 4096)
SCHEDULES = ("burst", "adversarial", "repair-only")
LOSSES = ("0.35", "0.50", "0.65")
CANDIDATES = tuple(range(256))
DEV_ROOTS = (
    "0xfdd68552844e91b2", "0x6163459a37898740",
    "0x3781d277770aeadc", "0xbcc1940397e96027",
    "0xa64ba14f42d6bd11", "0x066ca27133fd832a",
)
H1_ROOT_COUNT = 6
H2_ROOT_COUNT = 4
TIMING_CORE = 0
TIMING_SETUP_SCHEDULE = "burst"
DEV_SEED_LEDGER_SHA256 = (
    "6a6d6b979212a794258f7e37558c66ab61a35d7571d28015aceaa3b24f422580"
)

DRAND_CHAIN_HASH = (
    "52db9ba70e0cc0f6eaf7803dd07447a1f5477735fd3f661792ba94600c84e971"
)
DRAND_PUBLIC_KEY = (
    "83cf0f2896adee7eb8b5f01fcad3912212c437e0073e911fb90022d3e760183c8"
    "c4b450b6a0a6c3ac6a5776a2d1064510d1fec758c921cc22b0e17e63aaf4bcb5"
    "ed66304de9cf809bd274ca73bab4af5a6e9c76a4bc09e76eae8991ef5ece45a"
)
DRAND_GROUP_HASH = (
    "f477d5c89f21a17c863a7f937c6a6d15859414d2be09cd448d4279af331c5d3e"
)
DRAND_GENESIS_TIME = 1_692_803_367
DRAND_PERIOD_SECONDS = 3
DRAND_RELAYS = (
    "api.drand.sh", "api2.drand.sh", "api3.drand.sh",
    "drand.cloudflare.com",
)
DRAND_CLIENT_VERSION = "1.4.2"
DRAND_CLIENT_INTEGRITY = (
    "sha512-jeNJmrVplfgIA/GVndxxJ5mo8y63BS2pEdNhk1siU4pQ+z/BnxsqRnxjH9ag"
    "1ip887s12SEgo0MTZPbQNz27NA=="
)
DRAND_CLIENT_TARBALL_SHA1 = "f9108eef6881e62c0c0f154f30f7bd0a818ea809"
DRAND_CLIENT_TAG_COMMIT = "ef8c9260294f8699b5e8c27a6b764f8f0d768bea"
DRAND_CLIENT_TARBALL_SHA256 = (
    "81de34afba38520b461152bf032cfb5139bb6ced205bf9f50bc8216fdc394eef"
)
DRAND_CLIENT_TARBALL_BYTES = 551_599
DRAND_CLIENT_BUNDLE_MEMBER = "package/build/cjs/index.cjs"
DRAND_CLIENT_BUNDLE_SHA256 = (
    "45cb65d533cc7e8527e9bba92df875c066511c3d6286adc7fcb293f0d03c7566"
)
DRAND_CLIENT_BUNDLE_BYTES = 277_729
DRAND_NODE_VERSION = "v18.19.1"
DRAND_NPM_VERSION = "9.2.0"
DRAND_CANONICAL_KNOWN_SHA256 = (
    "b03d725b13e383cb643e02e448bf663e70e0f1b9dc85c33b0e34cefe151efeff"
)
DRAND_MASTER_ROOT_KNOWN_SHA256 = (
    "631985ad4fbcf2f98b5b725d9115458eb1e7d34d173f96f65847f5f4c08d0ce5"
)
DRAND_KNOWN_ROUND = {
    "round": 1,
    "randomness": (
        "1466a6cd24e327188770752f6134001c64d6efcc590ccc26b721611ad96f165a"
    ),
    "signature": (
        "b55e7cb2d5c613ee0b2e28d6750aabbb78c39dcc96bd9d38c2c2e12198df955"
        "71de8e8e402a0cc48871c7089a2b3af4b"
    ),
}


@dataclass(frozen=True)
class RoundSpec:
    name: str
    input_count: int
    retain_count: int
    root_indexes: tuple[int, ...]
    widths: tuple[int, ...]
    expected_candidate_cells: int
    expected_jobs: int


ROUNDS = (
    RoundSpec("r1", 256, 64, (0,), (64,), 1_476_864, 9_216),
    RoundSpec("r2", 64, 16, (1,), (64, 1280), 738_432, 2_304),
    RoundSpec("r3", 16, 8, (2, 3), WIDTHS, 738_432, 2_304),
    RoundSpec("r4", 8, 1, (4, 5), WIDTHS, 369_216, 2_304),
)

EMPTY_SHA256 = hashlib.sha256(b"").hexdigest()
MANIFEST_MAGIC = b"wirehair-wh2-holdout-seal-manifest-v1\n"


@dataclass(frozen=True)
class HoldoutPhaseSpec:
    phase: str
    table_name: str
    seal_tree: str
    root_file: str
    required_files: tuple[str, ...]
    excluded_prefixes: tuple[str, ...]
    forbidden_before_seal: tuple[str, ...]


HOLDOUT_PHASES = {
    "h1": HoldoutPhaseSpec(
        phase="h1",
        table_name="development",
        seal_tree="development_table",
        root_file="holdout/h1_roots.json",
        required_files=(
            "prepare.json", "frozen/contract.json", "frozen/protocol.json",
            "frozen/staged.sha256", "frozen/wirehair_v2_bench",
            "frozen/wh2_preferred_attempt_search.py",
            "frozen/wh2_preferred_attempt_campaign.py",
            "frozen/wh2_preferred_attempt_holdout.py",
            "frozen/wh2_preferred_attempt_timing.py",
            "frozen/wh2_drand_verify.cjs", "frozen/drand-client.cjs",
            "frozen/drand-package.json", "frozen/drand-package-lock.json",
            "frozen/r1_freeze_publication.json",
            "development/frozen_table.json",
            "development/phase_complete.sha256",
            "development/analysis_complete.sha256",
        ),
        excluded_prefixes=(
            "beacon/h1", "beacon/h2", "beacon/h1.lock", "beacon/h2.lock",
            "locks",
            "seals/development_table", "seals/h1_table",
            "holdout/h1_roots.json", "holdout/h2_roots.json", "h1", "h2",
            "timing", "thermal/h1-paired.csv", "thermal/h2-route.csv",
            "thermal/h2-paired.csv",
        ),
        forbidden_before_seal=(
            "holdout/h1_roots.json", "holdout/h2_roots.json", "h1", "h2",
            "timing", "thermal/h1-paired.csv", "thermal/h2-route.csv",
            "thermal/h2-paired.csv",
        ),
    ),
    "h2": HoldoutPhaseSpec(
        phase="h2",
        table_name="h1",
        seal_tree="h1_table",
        root_file="holdout/h2_roots.json",
        required_files=(
            "prepare.json", "frozen/contract.json", "frozen/protocol.json",
            "frozen/staged.sha256", "frozen/wirehair_v2_bench",
            "frozen/wh2_preferred_attempt_search.py",
            "frozen/wh2_preferred_attempt_campaign.py",
            "frozen/wh2_preferred_attempt_holdout.py",
            "frozen/wh2_preferred_attempt_timing.py",
            "frozen/wh2_drand_verify.cjs", "frozen/drand-client.cjs",
            "frozen/drand-package.json", "frozen/drand-package-lock.json",
            "frozen/r1_freeze_publication.json",
            "development/frozen_table.json",
            "development/phase_complete.sha256",
            "development/analysis_complete.sha256",
            "holdout/h1_roots.json", "h1/frozen_table.json",
            "h1/phase_complete.json", "h1/phase_complete.sha256",
            "h1/analysis_complete.json", "h1/analysis_complete.sha256",
            "h1/route_cache/index.json",
            "h1/route_cache/index.json.sha256",
            "h1/controller_complete.json",
            "h1/controller_complete.json.sha256",
        ),
        excluded_prefixes=(
            "beacon/h2", "beacon/h1.lock", "beacon/h2.lock", "locks",
            "seals/h1_table", "holdout/h2_roots.json", "h2", "timing",
            "thermal/h2-route.csv", "thermal/h2-paired.csv",
        ),
        forbidden_before_seal=(
            "holdout/h2_roots.json", "h2", "timing",
            "thermal/h2-route.csv", "thermal/h2-paired.csv",
        ),
    ),
}


class CampaignError(RuntimeError):
    pass


def die(message: str) -> None:
    raise CampaignError(message)


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ensure_ascii=True,
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def exact_decimal_milli(text: str, context: str) -> int:
    try:
        scaled = Decimal(text) * 1000
    except (InvalidOperation, ValueError):
        die(f"{context}: invalid decimal {text!r}")
    if (not scaled.is_finite() or scaled < 0 or
            scaled != scaled.to_integral_value()):
        die(f"{context}: value is not exact to 0.001")
    value = int(scaled)
    if value >= 1 << 64:
        die(f"{context}: value does not fit u64")
    return value


def strict_uint(text: str, context: str) -> int:
    if not text or (text != "0" and text.startswith("0")) or not text.isdigit():
        die(f"{context}: noncanonical unsigned integer {text!r}")
    return int(text, 10)


def development_seed_ledger() -> dict[str, list[str]]:
    """Return only roots that are allowed to exist before development."""
    ledger = {"dev": list(DEV_ROOTS)}
    encoded = (canonical_json(ledger) + "\n").encode("ascii")
    if sha256_bytes(encoded) != DEV_SEED_LEDGER_SHA256:
        die("development seed ledger hash changed")
    values = [int(seed, 0) for roots in ledger.values() for seed in roots]
    if len(values) != len(set(values)) or any(not 0 <= value < 1 << 64 for value in values):
        die("development seed ledger is not a distinct u64 domain")
    return ledger


def future_beacon_contract() -> dict[str, Any]:
    return {
        "provider": "League of Entropy drand Quicknet",
        "chain": {
            "hash": DRAND_CHAIN_HASH,
            "public_key": DRAND_PUBLIC_KEY,
            "group_hash": DRAND_GROUP_HASH,
            "scheme_id": "bls-unchained-g1-rfc9380",
            "beacon_id": "quicknet",
            "genesis_time_utc_s": DRAND_GENESIS_TIME,
            "period_seconds": DRAND_PERIOD_SECONDS,
        },
        "relays": [f"https://{host}" for host in DRAND_RELAYS],
        "api": f"/<chain-hash>/public/<exact-decimal-round>",
        "latest_endpoint_for_entropy": "forbidden",
        "verified_latest_preseal_lower_bound": {
            "role": "clock/liveness lower bound only; never root entropy",
            "endpoint": f"/<chain-hash>/public/latest",
            "verification": "same pinned four-origin BLS and >=3 quorum",
            "record": "latest round, canonical hash, raw response hashes, quorum",
            "failure": "no seal is created",
        },
        "selection": {
            "h1_preceding_seal": "development_table_seal.json",
            "h2_preceding_seal": "h1_table_seal.json",
            "threshold": "D=seal_ms+120000 using BigInt integer milliseconds",
            "round": (
                "R=max(1+ceil((D-G)/3000), verified_latest_R+40, "
                "abandoned_R+1); T=G+(R-1)*3000; assert T>=D, "
                "R>verified_latest_R, R>abandoned_R, and R is a JavaScript "
                "safe integer"
            ),
            "seal_publication_deadline": (
                "canonical seal record and manifest must be atomically fsynced "
                "and its hash published before T"
            ),
            "reseal": (
                "an unconfirmed or late attempt is permanently abandoned and "
                "a fresh seal/tag/round may be created only without querying "
                "the abandoned R"
            ),
            "unavailable_invalid_or_no_three_relay_consensus": (
                "block and retry the same exact round only"
            ),
            "fallback_substitution_or_later_round": "forbidden",
        },
        "verification": {
            "minimum_identical_relays": 3,
            "provider_diversity": (
                "three api.drand.sh relays plus drand.cloudflare.com"
            ),
            "transport": (
                "custom pinned ChainClient; no remote info trust; HTTPS GET; "
                "8 second timeout; redirects forbidden; 4096-byte response cap"
            ),
            "each_relay": (
                "exact round; pinned full chain info; randomness equals "
                "SHA256(signature); RFC9380 BLS signature under pinned key"
            ),
            "agreement": "round, randomness, and signature byte-identical",
            "verified_disagreement": "permanent impossible-condition block",
            "canonical_beacon": (
                "exact ASCII fixed-key-order {round,randomness,signature} plus LF; "
                "hex lowercased; unrelated response fields ignored"
            ),
            "canonical_round1_sha256": DRAND_CANONICAL_KNOWN_SHA256,
            "wrapper": "frozen/wh2_drand_verify.cjs",
            "known_round_self_test": DRAND_KNOWN_ROUND,
        },
        "verifier_supply_chain": {
            "package": f"drand-client@{DRAND_CLIENT_VERSION}",
            "registry": "https://registry.npmjs.org",
            "tarball": (
                "https://registry.npmjs.org/drand-client/-/"
                f"drand-client-{DRAND_CLIENT_VERSION}.tgz"
            ),
            "npm_integrity": DRAND_CLIENT_INTEGRITY,
            "tarball_sha1": DRAND_CLIENT_TARBALL_SHA1,
            "git_tag_commit": DRAND_CLIENT_TAG_COMMIT,
            "tarball_bytes": DRAND_CLIENT_TARBALL_BYTES,
            "tarball_sha256": DRAND_CLIENT_TARBALL_SHA256,
            "bundle_member": DRAND_CLIENT_BUNDLE_MEMBER,
            "bundle_bytes": DRAND_CLIENT_BUNDLE_BYTES,
            "bundle_sha256": DRAND_CLIENT_BUNDLE_SHA256,
            "node_version": DRAND_NODE_VERSION,
            "npm_version": DRAND_NPM_VERSION,
            "install": "npm ci --ignore-scripts --no-audit --no-fund",
            "package_scripts": "disabled",
            "frozen_artifacts": [
                "frozen/drand-package.json",
                "frozen/drand-package-lock.json",
                "frozen/drand-client.cjs",
            ],
        },
        "root_derivation": {
            "counts": {"h1": H1_ROOT_COUNT, "h2": H2_ROOT_COUNT},
            "master_input": (
                "ASCII wirehair-wh2-holdout-root-v1 + NUL, preceding-seal "
                "SHA256 raw32, chain hash raw32, R u64be, signature raw48, "
                "randomness raw32; exactly 181 bytes"
            ),
            "master_output": "SHA256 of master input",
            "round1_zero_seal_master_golden": DRAND_MASTER_ROOT_KNOWN_SHA256,
            "expansion": (
                "SHA256(ASCII phase-domain + NUL + master-root raw32 + index "
                "u32be), first 8 bytes big-endian u64"
            ),
            "phase_domains": {
                "h1": "wirehair-wh2-preferred-attempt-h1-seed-v1",
                "h2": "wirehair-wh2-preferred-attempt-h2-seed-v1",
            },
            "round1_zero_seal_expansion_golden": {
                "h1": [
                    "0xa16a872882abdc61", "0x04e0a1e6666926c6",
                    "0x9276294f95165cf2", "0xca3a62bcae7dff7d",
                    "0xa04c420c49e6946f", "0xc89a8bc53784d3f7",
                ],
                "h2": [
                    "0xf8be77afe18e6b43", "0xf4588ddd7097dfae",
                    "0xec4ea66c4756638b", "0x298cf149944eeafa",
                ],
            },
            "distinctness": "within phase and from every earlier root",
        },
        "controller": {
            "states": [
                "SEALED", "WAITING_UNTIL_T", "WAITING_FOR_QUORUM",
                "ROOTED", "ABANDONED_*", "BLOCKED_PERMANENT",
            ],
            "wave_origins": 4,
            "retry_delays_seconds": [1, 2, 4, 8, 16, 30],
            "retry_tail_seconds": 30,
            "temporary_exit_after_seconds": 900,
            "temporary_exit": "persist WAITING and resume the exact same seal/R",
            "permanent_blocks": [
                "sealed file, lock, helper, executable, or runtime hash mutation",
                "target formula inconsistency or clock rollback before seal",
                "BLS-verified disagreement", "stored-root mismatch",
            ],
            "phase_manifest": {
                "magic": MANIFEST_MAGIC.decode("ascii").rstrip("\n"),
                "record": "sha256lower SP decimal-byte-size SP relative-POSIX-path LF",
                "order": "raw normalized relative path byte order",
                "validation": (
                    "NFC; no absolute/.., CR/LF/NUL/backslash, duplicate, "
                    "symlink, nonregular, missing required, or future-phase file"
                ),
                "h1_required": list(HOLDOUT_PHASES["h1"].required_files),
                "h2_required": list(HOLDOUT_PHASES["h2"].required_files),
            },
            "publication": {
                "remote": "origin on github.com",
                "tag": "wh2-h12-kpreferred-{phase}-seal-{seal_sha256[0:16]}",
                "kind": "unique annotated non-force tag pointing to source commit",
                "annotation_binds": [
                    "full seal-record SHA256", "manifest SHA256", "phase",
                    "target R", "target T", "source commit",
                ],
                "deadline": "push and exact ls-remote tag+peeled verification before T",
                "late_or_failed": (
                    "retain ABANDONED audit/tag; never query, delete, force, or "
                    "reuse; create a fresh seal/tag/round"
                ),
            },
            "r1_freeze_publication": {
                "deadline": "confirmed before any R1 recovery observation",
                "tag": "wh2-h12-kpreferred-r1-freeze-{prepare_sha256[0:16]}",
                "kind": "unique annotated non-force tag at source commit",
                "binds": [
                    "prepare SHA256", "staged-seal SHA256", "source commit",
                    "binary/controller/contract/protocol/route-context SHA256",
                    "development seed ledger and future-beacon contract SHA256",
                ],
                "mutation": "requires a new prepare record and new tag",
            },
        },
        "root_file_state_machine": {
            "before_development_seal": [
                "holdout/h1_roots.json absent", "holdout/h2_roots.json absent",
            ],
            "after_development_seal": (
                "only h1_roots.json may be acquired from its exact target pulse"
            ),
            "before_h1_table_seal": "holdout/h2_roots.json absent",
            "after_h1_table_seal": (
                "only h2_roots.json may be acquired from its exact target pulse"
            ),
        },
    }


def target_beacon_round(seal_ms: int) -> int:
    if not isinstance(seal_ms, int) or isinstance(seal_ms, bool) or seal_ms < 0:
        die("seal UTC timestamp must be a nonnegative integer millisecond count")
    threshold_ms = seal_ms + 120_000
    genesis_ms = DRAND_GENESIS_TIME * 1000
    period_ms = DRAND_PERIOD_SECONDS * 1000
    if threshold_ms < genesis_ms:
        die("seal threshold predates the pinned Quicknet genesis")
    intervals = (threshold_ms - genesis_ms + period_ms - 1) // period_ms
    round_number = intervals + 1
    scheduled_ms = genesis_ms + (round_number - 1) * period_ms
    if (not genesis_ms <= threshold_ms <= scheduled_ms < threshold_ms + period_ms or
            round_number > (1 << 53) - 1):
        die("internal Quicknet target-round arithmetic mismatch")
    return round_number


def canonical_beacon_bytes(beacon: dict[str, Any]) -> bytes:
    if not isinstance(beacon, dict):
        die("verified beacon must be an object")
    round_number = beacon.get("round")
    signature = beacon.get("signature")
    randomness = beacon.get("randomness")
    if (set(beacon) != {"round", "signature", "randomness"} or
            not isinstance(round_number, int) or isinstance(round_number, bool) or
            not 1 <= round_number < 1 << 64 or
            not isinstance(signature, str) or
            not re.fullmatch(r"[0-9a-f]{96}", signature) or
            not isinstance(randomness, str) or
            not re.fullmatch(r"[0-9a-f]{64}", randomness) or
            hashlib.sha256(bytes.fromhex(signature)).hexdigest() != randomness):
        die("verified beacon fields are not canonical or self-consistent")
    return (
        f'{{"round":{round_number},"randomness":"{randomness}",'
        f'"signature":"{signature}"}}\n'
    ).encode("ascii")


def holdout_master_root(
    preceding_seal_sha256: str,
    beacon: dict[str, Any],
) -> str:
    if not re.fullmatch(r"[0-9a-f]{64}", preceding_seal_sha256):
        die("preceding seal SHA256 is not canonical")
    canonical = canonical_beacon_bytes(beacon)
    round_number = beacon["round"]
    signature = beacon["signature"]
    randomness = beacon["randomness"]
    if (round_number == 1 and signature == DRAND_KNOWN_ROUND["signature"] and
            randomness == DRAND_KNOWN_ROUND["randomness"] and
            sha256_bytes(canonical) != DRAND_CANONICAL_KNOWN_SHA256):
        die("canonical Quicknet round-one golden changed")
    master_input = (
        b"wirehair-wh2-holdout-root-v1\0" +
        bytes.fromhex(preceding_seal_sha256) + bytes.fromhex(DRAND_CHAIN_HASH) +
        round_number.to_bytes(8, "big") + bytes.fromhex(signature) +
        bytes.fromhex(randomness)
    )
    if len(master_input) != 181:
        die("holdout master-root input length changed")
    master_root = hashlib.sha256(master_input).digest()
    if (preceding_seal_sha256 == "0" * 64 and round_number == 1 and
            signature == DRAND_KNOWN_ROUND["signature"] and
            randomness == DRAND_KNOWN_ROUND["randomness"] and
            master_root.hex() != DRAND_MASTER_ROOT_KNOWN_SHA256):
        die("holdout master-root golden changed")
    return master_root.hex()


def derive_beacon_roots(
    phase: str,
    preceding_seal_sha256: str,
    beacon: dict[str, Any],
    earlier_holdout_roots: Sequence[str] = (),
) -> tuple[str, ...]:
    counts = {"h1": H1_ROOT_COUNT, "h2": H2_ROOT_COUNT}
    if phase not in counts:
        die(f"invalid holdout phase {phase!r}")
    master_root = bytes.fromhex(
        holdout_master_root(preceding_seal_sha256, beacon))
    phase_domain = f"wirehair-wh2-preferred-attempt-{phase}-seed-v1".encode("ascii")
    roots = tuple(
        "0x" + hashlib.sha256(
            phase_domain + b"\0" + master_root + index.to_bytes(4, "big")
        ).digest()[:8].hex()
        for index in range(counts[phase])
    )
    if ((phase == "h1" and earlier_holdout_roots) or
            (phase == "h2" and len(earlier_holdout_roots) != H1_ROOT_COUNT)):
        die("earlier holdout root coverage is invalid for this phase")
    try:
        earlier = {
            int(root, 0) for root in (*DEV_ROOTS, *earlier_holdout_roots)
            if re.fullmatch(r"0x[0-9a-f]{16}", root)
        }
    except (TypeError, ValueError):
        die("earlier holdout roots are not canonical u64 values")
    if len(earlier) != len(DEV_ROOTS) + len(earlier_holdout_roots):
        die("earlier holdout roots are duplicate or noncanonical")
    values = [int(root, 0) for root in roots]
    if len(set(values)) != len(values) or earlier.intersection(values):
        die("derived holdout roots collide with an earlier root")
    return roots


def validate_beacon_contract_kats() -> None:
    genesis_ms = DRAND_GENESIS_TIME * 1000
    for delta, expected in ((0, 1), (1, 2), (2999, 2), (3000, 2), (3001, 3)):
        if target_beacon_round(genesis_ms + delta - 120_000) != expected:
            die(f"Quicknet target-round boundary KAT failed at delta {delta}")
    zero_seal = "0" * 64
    h1_expected = (
        "0xa16a872882abdc61", "0x04e0a1e6666926c6",
        "0x9276294f95165cf2", "0xca3a62bcae7dff7d",
        "0xa04c420c49e6946f", "0xc89a8bc53784d3f7",
    )
    h2_expected = (
        "0xf8be77afe18e6b43", "0xf4588ddd7097dfae",
        "0xec4ea66c4756638b", "0x298cf149944eeafa",
    )
    h1 = derive_beacon_roots("h1", zero_seal, DRAND_KNOWN_ROUND)
    h2 = derive_beacon_roots(
        "h2", zero_seal, DRAND_KNOWN_ROUND, earlier_holdout_roots=h1)
    if h1 != h1_expected or h2 != h2_expected:
        die("holdout phase-expansion golden changed")


def verify_source_seals(source_root: Path) -> None:
    expected = {
        source_root / "analysis/failures_and_shortfalls.csv":
            SOURCE_FAILURES_SHA256,
        source_root / "analysis/paired_cells.csv": SOURCE_PAIRED_SHA256,
        source_root / "results/phase_complete.sha256":
            SOURCE_PHASE_SEAL_SHA256,
        source_root / "analysis/analysis_complete.sha256":
            SOURCE_ANALYSIS_SEAL_SHA256,
    }
    for path, digest in expected.items():
        if common.sha256_file(path.resolve(strict=True)) != digest:
            die(f"sealed source hash mismatch: {path}")
    source_contract_path = source_root / "frozen/contract.json"
    if common.sha256_file(source_contract_path) != SOURCE_CONTRACT_SHA256:
        die("sealed source campaign contract digest mismatch")
    contract = json.loads(common.stable_bytes(source_contract_path))
    if (not isinstance(contract, dict) or
            contract.get("source_commit") != SOURCE_COMMIT or
            contract.get("block_bytes") != 64 or
            contract.get("loss") != "0.50" or
            contract.get("schedules") != list(SCHEDULES) or
            contract.get("packet_peel_seed_xor") != "0x0" or
            contract.get("trials") != 1):
        die("sealed source campaign contract mismatch")
    # Verify every source result and every published analysis file before using
    # either outcome or timing data.  The hashes above pin the manifests too.
    common.verify_sha_manifest(
        source_root, source_root / "results/phase_complete.sha256")
    common.verify_sha_manifest(
        source_root, source_root / "analysis/analysis_complete.sha256")


def derive_active_cohort(source_root: Path) -> tuple[int, ...]:
    path = source_root / "analysis/failures_and_shortfalls.csv"
    required = {
        "K", "active_packet_peel_seed_xor",
        "two_anchor_adaptive_rank_fail", "two_anchor_adaptive_error",
    }
    seen_failure_cells: set[tuple[int, str, str]] = set()
    active: set[int] = set()
    with io.StringIO(common.stable_bytes(path).decode("utf-8"), newline="") as text:
        reader = csv.DictReader(text)
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            die("sealed failure source has an unexpected schema")
        for row in reader:
            K = strict_uint(row["K"], "source K")
            rank_fail = strict_uint(
                row["two_anchor_adaptive_rank_fail"], "source rank_fail")
            error = strict_uint(
                row["two_anchor_adaptive_error"], "source error")
            if not K_MIN <= K <= K_MAX or rank_fail > 1 or error > 1:
                die("sealed failure source has an invalid outcome")
            if error:
                die("sealed q0 cohort source contains a codec error")
            if K < CUTOFF or not rank_fail:
                continue
            if int(row["active_packet_peel_seed_xor"], 0) != 0:
                die("active q0 cohort contains a nonzero packet-seed xor")
            key = (K, row.get("seed", ""), row.get("schedule", ""))
            if key in seen_failure_cells:
                die(f"duplicate q0 failure cell {key}")
            seen_failure_cells.add(key)
            active.add(K)
    cohort = tuple(sorted(active))
    encoded = "".join(f"{K}\n" for K in cohort).encode("ascii")
    multiplicity: dict[int, int] = {K: 0 for K in cohort}
    for K, _seed, _schedule in seen_failure_cells:
        multiplicity[K] += 1
    if (len(cohort) != ACTIVE_K_COUNT or len(seen_failure_cells) != 656 or
            cohort[0] != 4272 or cohort[-1] != 63816 or
            sum(count > 1 for count in multiplicity.values()) != 14 or
            max(multiplicity.values()) != 3 or len(encoded) != 3823 or
            sha256_bytes(encoded) != COHORT_SHA256):
        die("derived active q0 cohort contract changed")
    return cohort


def source_timing_weights(
    source_root: Path, cohort: Sequence[int],
) -> tuple[dict[int, int], dict[int, int]]:
    cohort_set = set(cohort)
    active_weights = {K: 0 for K in cohort}
    active_counts = {K: 0 for K in cohort}
    all_weights = {K: 0 for K in range(K_MIN, K_MAX + 1)}
    all_counts = {K: 0 for K in range(K_MIN, K_MAX + 1)}
    seen_cells: set[tuple[int, int, str]] = set()
    stdout_root = source_root / "results/stdout"
    paths = sorted(stdout_root.glob("*.two_anchor_adaptive.*.csv"))
    if len(paths) != 2781:
        die(f"sealed q0 timing source has {len(paths)} files, want 2781")
    for path in paths:
        match = re.fullmatch(
            r"job[0-9]{5}\.two_anchor_adaptive\.seed([0-2])\."
            r"(burst|adversarial|repair-only)\.group[0-9]{3}\."
            r"(small|large)\.csv",
            path.name,
        )
        if match is None:
            die(f"q0 timing filename mismatch: {path.name}")
        seed_index = int(match.group(1))
        schedule = match.group(2)
        size_class = match.group(3)
        data = common.stable_bytes(path).decode("utf-8")
        first, separator, csv_text = data.partition("\n")
        tokens = first.split()
        if not separator or tokens[:2] != ["#", "precodefail:"]:
            die(f"q0 timing preamble mismatch: {path.name}")
        preamble: dict[str, str] = {}
        for token in tokens[2:]:
            if token.count("=") != 1:
                die(f"q0 timing preamble token mismatch: {path.name}")
            key, value = token.split("=", 1)
            if not key or not value or key in preamble:
                die(f"q0 timing duplicate/empty preamble token: {path.name}")
            preamble[key] = value
        expected_preamble = {
            "trials": "1", "loss": "0.5",
            "seed": SOURCE_ROOTS[seed_index], "completion": "mixed",
            "mixed_period": "244", "mixed_gf256_rows": "10",
            "mixed_gf16_rows": "2", "mixed_geometry": "frozen",
            "mixed_residue_schedule": "constant",
            "packet_peel_seed_xor": "0x0",
            "binary_dense_two_anchor": "0" if size_class == "small" else "1",
            "schedule": schedule,
        }
        if (any(preamble.get(key) != value
                for key, value in expected_preamble.items()) or
                "mix_count" in preamble):
            die(f"q0 timing preamble mismatch: {path.name}")
        reader = csv.DictReader(io.StringIO(csv_text, newline=""))
        required = {
            "N", "bb", "mix_count", "overhead", "trials",
            "solve_ms_mu", "active_packet_peel_seed_xor",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            die(f"q0 timing CSV schema mismatch: {path.name}")
        rows = 0
        for row in reader:
            rows += 1
            K = strict_uint(row["N"], "timing K")
            cell_key = (K, seed_index, schedule)
            if cell_key in seen_cells:
                die(f"duplicate q0 timing cell: {cell_key}")
            seen_cells.add(cell_key)
            if (not K_MIN <= K <= K_MAX or row["bb"] != "64" or
                    row["mix_count"] != "2" or row["overhead"] != "0" or
                    row["trials"] != "1" or
                    int(row["active_packet_peel_seed_xor"], 0) != 0 or
                    (K < CUTOFF) != (size_class == "small")):
                die(f"q0 timing row geometry mismatch: {path.name}/{K}")
            micros = exact_decimal_milli(
                row["solve_ms_mu"], f"{path.name}/{K}/solve_ms_mu")
            all_weights[K] += micros
            all_counts[K] += 1
            if K in cohort_set:
                active_weights[K] += micros
                active_counts[K] += 1
        if rows == 0:
            die(f"empty q0 timing source: {path.name}")
    if set(all_counts.values()) != {9} or set(active_counts.values()) != {9}:
        die("q0 timing source does not cover exactly nine cells per K")
    expected_cells = {
        (K, seed_index, schedule)
        for K in range(K_MIN, K_MAX + 1)
        for seed_index in range(len(SOURCE_ROOTS))
        for schedule in SCHEDULES
    }
    if seen_cells != expected_cells:
        die("q0 timing source does not cover the exact K/root/schedule cube")
    return active_weights, all_weights


def lpt_bins(
    domain: Sequence[int], weights: dict[int, int],
) -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...], dict[int, int]]:
    heap: list[tuple[int, int]] = [(0, index) for index in range(BIN_COUNT)]
    heapq.heapify(heap)
    contents: list[list[int]] = [[] for _ in range(BIN_COUNT)]
    assignment: dict[int, int] = {}
    for K in sorted(domain, key=lambda value: (-weights[value], value)):
        load, index = heapq.heappop(heap)
        contents[index].append(K)
        assignment[K] = index
        heapq.heappush(heap, (load + weights[K], index))
    loads = [0] * BIN_COUNT
    for K, index in assignment.items():
        loads[index] += weights[K]
    return (
        tuple(tuple(sorted(values)) for values in contents),
        tuple(loads), assignment,
    )


def serialize_source_ledgers(
    cohort: Sequence[int], active_weights: dict[int, int],
    all_weights: dict[int, int],
) -> dict[str, bytes]:
    bins, loads, assignment = lpt_bins(cohort, active_weights)
    all_domain = tuple(range(K_MIN, K_MAX + 1))
    all_bins, all_loads, all_assignment = lpt_bins(all_domain, all_weights)
    ledgers = {
        "cohort.txt": "".join(f"{K}\n" for K in cohort).encode("ascii"),
        "weights.tsv": "".join(
            f"{K}\t{active_weights[K]}\t9\n" for K in cohort
        ).encode("ascii"),
        "assignment.tsv": "".join(
            f"{K}\t{assignment[K]}\t{active_weights[K]}\n" for K in cohort
        ).encode("ascii"),
        "bins.tsv": "".join(
            f"{index}\t{loads[index]}\t{len(bins[index])}\t"
            f"{','.join(map(str, bins[index]))}\n"
            for index in range(BIN_COUNT)
        ).encode("ascii"),
        "allk_assignment.tsv": "".join(
            f"{K}\t{all_assignment[K]}\t{all_weights[K]}\n"
            for K in all_domain
        ).encode("ascii"),
        "allk_bins.tsv": "".join(
            f"{index}\t{all_loads[index]}\t{len(all_bins[index])}\t"
            f"{','.join(map(str, all_bins[index]))}\n"
            for index in range(BIN_COUNT)
        ).encode("ascii"),
    }
    expected = {
        "cohort.txt": (3823, COHORT_SHA256),
        "weights.tsv": (9786, WEIGHTS_SHA256),
        "assignment.tsv": (10517, ASSIGNMENT_SHA256),
        "bins.tsv": (5505, BINS_SHA256),
        "allk_assignment.tsv": (1017797, ALLK_ASSIGNMENT_SHA256),
        "allk_bins.tsv": (375086, ALLK_BINS_SHA256),
    }
    for name, (size, digest) in expected.items():
        if len(ledgers[name]) != size or sha256_bytes(ledgers[name]) != digest:
            die(f"derived {name} fingerprint mismatch")
    if (sum(loads) != 511_795_298 or min(loads) != 3_967_753 or
            max(loads) != 4_055_182 or
            min(map(len, bins)) != 4 or max(map(len, bins)) != 6 or
            sum(all_loads) != 34_236_201_584 or
            min(all_loads) != 267_469_573 or
            max(all_loads) != 267_470_677 or
            min(map(len, all_bins)) != 499 or max(map(len, all_bins)) != 500):
        die("derived LPT bin arithmetic changed")
    return ledgers


def protocol_definition() -> dict[str, Any]:
    validate_beacon_contract_kats()
    seeds = development_seed_ledger()
    beacon = future_beacon_contract()
    round_records = []
    for spec in ROUNDS:
        calculated = (
            ACTIVE_K_COUNT * spec.input_count * len(spec.root_indexes) *
            len(SCHEDULES) * len(LOSSES) * len(spec.widths)
        )
        if calculated != spec.expected_candidate_cells:
            die(f"{spec.name} candidate arithmetic mismatch")
        batches = (spec.input_count + 31) // 32
        jobs = (
            BIN_COUNT * batches * len(spec.root_indexes) *
            len(SCHEDULES) * len(LOSSES)
        )
        if jobs != spec.expected_jobs:
            die(f"{spec.name} job arithmetic mismatch")
        route_manifest_jobs = BIN_COUNT * batches
        round_records.append({
            "name": spec.name,
            "input_candidates_per_K": spec.input_count,
            "retain_candidates_per_K": spec.retain_count,
            "dev_root_indexes": list(spec.root_indexes),
            "widths": list(spec.widths),
            "schedule_loss_strata": 9,
            "candidate_cells": calculated,
            "candidate_batches_per_K": batches,
            "jobs": jobs,
            "route_manifest_jobs_max": route_manifest_jobs,
        })
    candidate_cells = sum(row["candidate_cells"] for row in round_records)
    candidate_jobs = sum(row["jobs"] for row in round_records)
    control_cells = ACTIVE_K_COUNT * len(DEV_ROOTS) * 9 * len(WIDTHS)
    control_jobs = BIN_COUNT * len(DEV_ROOTS) * 9
    h1_logical_cells = 2 * ACTIVE_K_COUNT * H1_ROOT_COUNT * 9 * len(WIDTHS)
    h1_jobs = BIN_COUNT * H1_ROOT_COUNT * 9
    h2_logical_cells = 2 * (K_MAX - K_MIN + 1) * H2_ROOT_COUNT * 9
    h2_jobs = BIN_COUNT * H2_ROOT_COUNT * 9
    development_route_jobs = sum(
        row["route_manifest_jobs_max"] for row in round_records)
    h2_route_jobs = BIN_COUNT
    if (
        candidate_cells != 3_322_944 or control_cells != 138_456 or
        candidate_cells + control_cells != 3_461_400 or
        candidate_jobs != 16_128 or control_jobs != 6_912 or
        h1_logical_cells != 276_912 or h1_jobs != 6_912 or
        h2_logical_cells != 4_607_928 or h2_jobs != 4_608 or
        candidate_jobs + control_jobs + h1_jobs + h2_jobs != 34_560
    ):
        die("global protocol arithmetic mismatch")
    return {
        "schema": "wirehair.wh2.h12_preferred_attempt.protocol.v1",
        "purpose": "post-selection H12 q0 seed-route tuning only",
        "production_changes_before_frozen_winners": False,
        "profile_id_before_frozen_winners": False,
        "generator_contracts": {"precode": 4, "packet_rows": 4},
        "future_seed_route_contract": {
            "policy": "KPreferredAttemptV1", "version": 1,
            "status": "forbidden until H2 and timing acceptance",
        },
        "K_domain": [K_MIN, K_MAX],
        "active_K": ACTIVE_K_COUNT,
        "active_K_predicate": (
            "sealed two_anchor_adaptive K>=4096 xor=0; errors must be zero; "
            "any rank_fail over nine bb64/loss.50/trials1 cells"
        ),
        "block_widths": list(WIDTHS),
        "schedules": list(SCHEDULES),
        "losses": list(LOSSES),
        "candidate_attempts": [0, 255],
        "candidate_batch": 32,
        "bins": BIN_COUNT,
        "development_roots": seeds,
        "holdout_roots_present": False,
        "future_holdout_beacon": beacon,
        "external_seed_derivation": {
            "dev": "sha256(wirehair.wh2.h12-preferred-attempt.dev-root.v1/i)[0:8]",
            "byte_order": "big-endian display; parsed as u64",
            "development_ledger_sha256": DEV_SEED_LEDGER_SHA256,
            "h1_h2": "see future_holdout_beacon.root_derivation",
        },
        "routing": {
            "canonical": "ascending first systematic-valid attempt a0 in 0..255",
            "preferred": "try p once at actual K and BlockBytes",
            "invalid": "fallback to canonical ascending selector skipping p",
            "direct": "preferred_valid && p != a0",
            "alias": "invalid fallback or p==a0; exact cached control observation",
            "logical_vs_physical": "aliases are logical cells and zero candidate recovery solves",
            "hard_eligibility": [
                "direct at bb64", "no fallback at any encountered width through H1",
                "direct at one or more encountered widths",
            ],
            "noncanonical_noop": "neutral alias, never an independent win",
            "manifest": {
                "schema": "preferredattempt-route-manifest.v1",
                "keys": "exact requested K x approved-width cross-product",
                "route_status": ["preferred", "control"],
                "control_status": (
                    "mandatory for every nonrouted H2 K; tuned arm is an exact "
                    "zero-recovery-solve alias of the one physical control solve"
                ),
                "context_binding": (
                    "SHA256 of canonical route_context.json containing source "
                    "commit, frozen binary/protocol/cohort/bin/root hashes"
                ),
                "probe_accounting": "additive exactly once per physical solve",
                "line_endings": "canonical LF (binary stdout on Windows)",
            },
        },
        "work_metric_symbols": {
            "X": "block_xors",
            "M": "block_muladds",
            "A": "inactivated columns",
        },
        "development": {
            "arithmetic_kind": "precommitted maxima; realized work may only shrink",
            "rounds": round_records,
            "candidate_cells": candidate_cells,
            "full_control_cube_cells": control_cells,
            "logical_cells": candidate_cells + control_cells,
            "candidate_jobs": candidate_jobs,
            "control_jobs": control_jobs,
            "jobs": candidate_jobs + control_jobs,
            "future_control_masking": (
                "a round may consume only its declared root/width panels even "
                "though d0..d5 controls are physically cached first"
            ),
            "survivor_ledger": (
                "per-K ordered eligible p IDs; keep up to 64/16/8/1 without "
                "backfill; a later-invalid p is removed, realized work shrinks, "
                "and a K with no survivor routes to control; canonical JSON and "
                "SHA256 sealed before the next panel"
            ),
            "hard_filter_before_ranking": [
                "any invariant or codec error",
                "not direct at bb64",
                "fallback at any width encountered through the current round",
            ],
            "per_K_lexicographic_rank": [
                "introductions_vs_matched_a0",
                "count adverse (bb,schedule,loss) aggregate slices vs matched a0",
                "max positive aggregate-slice failure delta",
                "sum positive aggregate-slice failure deltas",
                "total_candidate_failures",
                "negative_repairs_vs_matched_a0",
                "heavy_shortfall_sum",
                "binary_deficit_max",
                "binary_deficit_sum",
                "common-success_block_muladds_sum",
                "common-success_block_xors_sum",
                "inactivation_sum",
                "preferred_attempt_id",
            ],
            "ranking_panel_scope": (
                "cumulative over every declared matched panel through the "
                "current round; a repair/introduction is one exact "
                "candidate-vs-control logical cell"
            ),
            "r4_route_rule": (
                "route the sole survivor only when its cumulative rank tuple "
                "is strictly lexicographically better than the matched "
                "control tuple; ties route to control"
            ),
            "global_table_objective_report": [
                "errors", "max_K_failure_multiplicity", "multi_failure_K_count",
                "total_failures", "introductions", "worst_adverse_axis_delta",
                "sum_positive_axis_deltas", "heavy_shortfall_sum",
                "binary_deficit_gt15_count", "binary_deficit_max",
                "binary_deficit_sum", "common_success_xors",
                "common_success_muladds",
            ],
        },
        "h1": {
            "role": "sealed post-development pruning; never runner-up substitution",
            "roots": H1_ROOT_COUNT,
            "root_values_before_development_seal": "absent",
            "logical_cells": h1_logical_cells, "jobs": h1_jobs,
            "cells_per_K_per_arm": 216,
            "per_K_gates": [
                "seals valid; errors zero; direct at bb64; no fallback at four widths",
                "keep if repairs>introductions",
                "or repairs=introductions=0 and X/M/A all nonworse with at least one strict",
                "reject repairs=introductions>0",
                "each of three schedule totals and each of three loss totals nonworse",
                "36 width-by-schedule-by-loss slices report-only",
                "heavy shortfall sum, deficit>15 count, and max deficit nonworse",
                "for X/M/A: 200*candidate <= 201*control",
            ],
            "work_sum_domain": (
                "X/M/A integer sums use common-success matched cells only"
            ),
            "failure_action": "route K to canonical control; do not inspect runner-up",
        },
        "h2": {
            "role": "true final third-domain whole-table holdout; no edits",
            "roots": H2_ROOT_COUNT,
            "root_values_before_h1_table_seal": "absent",
            "width": 64,
            "logical_cells": h2_logical_cells, "jobs": h2_jobs,
            "gates": [
                "all seals valid and errors zero",
                "table failures strictly lower than control and repairs>introductions",
                "exact sign: 100*sum(comb(n,j),j=repairs..n) <= 2**n, n=repairs+introductions",
                "report nine schedule-loss cells; gate three schedule and three loss totals only",
                "weak-K, multi-K, max multiplicity, heavy shortfall, deficit>15, max deficit nonworse",
                "routed-cohort and all-K X/M/A: 200*candidate <= 201*control",
            ],
            "failure_multiplicity": (
                "weak K means >=1 rank_fail at bb64 over phase roots x nine "
                "schedule/loss strata; multi K means >=2; X/M/A integer sums "
                "use common-success matched cells only"
            ),
            "failure_action": "reject whole table; no edits or substitutions",
        },
        "timing": timing_contract(),
        "job_totals": {
            "development": candidate_jobs + control_jobs,
            "h1": h1_jobs, "h2": h2_jobs,
            "recovery_jobs_all_stages": 34_560,
            "route_manifest_jobs_max": development_route_jobs + h2_route_jobs,
            "all_process_jobs_max": (
                34_560 + development_route_jobs + h2_route_jobs
            ),
        },
        "thermal_edac_cleanup": {
            "workers": 128, "cpu_set": [0, 127], "balanced_bins": 128,
            "temperature_limits_c": {"cpu": 90.0, "dimm": 90.0},
            "timing_temperature_limits_c": {"cpu": 85.0, "dimm": 90.0},
            "consecutive_samples": 3,
            "stale_seconds": 5.0, "min_cpu_busy_pct": 95.0,
            "dimm_read_errors": 0, "edac_ce_delta": 0, "edac_ue_delta": 0,
            "child_process_groups_reaped": True,
        },
    }


def timing_contract() -> dict[str, Any]:
    return {
        "status": (
            "must pass after H1 table freeze and after the independent H2 "
            "root is irrevocably committed/revealed, before production promotion"
        ),
        "sample_priority": (
            "sha256('wirehair.wh2.h12-preferred-attempt.timing-sample.v1|' + decimal_K)"
        ),
        "bands": [[4096, 16383], [16384, 32767], [32768, 49151], [49152, 64000]],
        "selection": (
            "S=min(32,routed_K_count); up to 8 noncontrol K per band, then "
            "global hash fill; compare full SHA256 priority, ties by ascending "
            "numeric K; no substitutions; S<5 rejects, S=5..31 is valid"
        ),
        "minimum_for_exact_power": 5,
        "execution_core": TIMING_CORE,
        "trace_seed": {
            "domain": "wirehair.wh2.h12-preferred-attempt.timing-trace-seed.v1",
            "input": "domain + NUL + metric + NUL + decimal_K + NUL + decimal_bb + NUL + schedule",
            "output": "first 8 SHA256 bytes interpreted big-endian u64",
        },
        "solve_panels": {
            "widths": [64, 1280, 4096],
            "schedules": list(SCHEDULES), "loss": "0.50", "overhead": 4,
            "rhs": "distinct aligned prefaulted zero blocks at actual BlockBytes",
            "same_trace": True, "systems_and_runtimes_prebuilt_outside_timer": True,
            "cycle": ["C", "T", "T", "C", "T", "C", "C", "T"],
            "discarded_cycles": 1, "timed_cycles": 3,
            "runs_per_arm_panel": 16, "timed_runs_per_arm_panel": 12,
            "S32_total_solve_invocations": 9216,
            "S32_timed_solve_invocations": 6912,
            "realized_invocations": "total=288*S; timed=216*S",
        },
        "decoder_setup": {
            "operation": "BuildPrecodeSystem + PacketRowRuntime from transmitted actual attempt",
            "schedule": TIMING_SETUP_SCHEDULE,
            "excludes": "seed selection",
            "width_panels": [64, 1280, 4096],
            "same_cycle": True,
            "panel": "exact upper median of three timed cycle ratios",
            "K_ratio": "exact upper median of the three width-panel ratios",
            "global": "exact upper median of per-K setup ratios <=1.005",
            "each_width": (
                "exact upper median across sampled K <=1.005; "
                "200*numerator<=201*denominator"
            ),
            "S32_total_invocations": 3072,
            "S32_timed_invocations": 2304,
            "realized_invocations": "total=96*S; timed=72*S",
        },
        "statistics": {
            "upper_median": "sort ascending and select zero-based index floor(n/2)",
            "cycle_ratio": "sum(four T ns)/sum(four C ns)",
            "panel": "exact upper median of three cycle ratios",
            "K_solve": "exact upper median of nine panels",
            "global": (
                "upper median K ratio < 1 and exact one-sided sign p<=0.05"
            ),
            "sign_observations": (
                "per-K K_solve ratios; win iff T<C, loss iff T>C, exact ties "
                "excluded; n=wins+losses and n must be positive"
            ),
            "timing_ties": (
                "exact integer-nanosecond equality only; excluded from sign "
                "power but retained in medians and fully hashed ledgers"
            ),
            "exact_sign": (
                "20*sum(comb(n,j),j=wins..n) <= 2**n with integer arithmetic"
            ),
            "each_width": (
                "for each K take upper median of its three schedule panels, "
                "then upper median across K <=1.005; "
                "200*numerator<=201*denominator"
            ),
            "decoder_setup": "defined independently in decoder_setup",
        },
        "isolation": {
            "exclusive_physical_core": True, "same_numa_node": True,
            "smt_sibling_idle": True, "load_workers_stopped": True,
            "eviction_before_run": "max(2*LLC,256MiB)",
            "cycle_retry_limit": 2,
            "cycle_attempt_limit": 3,
            "cycle_retry_semantics": (
                "one original plus at most two replacement attempts for the "
                "exact contaminated cycle/panel; only predeclared contamination "
                "permits retry; never retry a clean result"
            ),
            "cycle_retry_record": (
                "each accepted or rejected cycle records SHA256 bindings for "
                "trace, arm order, setup/solve counters, raw nanoseconds, "
                "thermal interval, migration/fault/throttle evidence; retry "
                "the whole paired cycle only, never one arm selectively"
            ),
            "temperature_limits_c": {"cpu": 85.0, "dimm": 90.0},
            "contamination": [
                "migration", "major fault", "throttle", "APERF/MPERF excursion",
                "temperature excursion", "telemetry gap",
            ],
            "performance_monitor": {
                "source": "turbostat per-core APERF/MPERF-derived samples",
                "interval_ms": 100,
                "telemetry_gap": (
                    "no complete target-core sample, nonpositive Bzy/TSC MHz, "
                    "or TSC sample spread exceeds 1 MHz"
                ),
                "aperf_mperf_excursion": "100*max(Bzy_MHz) > 102*min(Bzy_MHz)",
                "throttle": "any 100*Bzy_MHz < 95*TSC_MHz",
                "raw_evidence_sha256_required": True,
            },
        },
        "panel_execution_order": (
            "K ascending; for each K solve widths 64,1280,4096 and schedules "
            "burst,adversarial,repair-only; then setup widths 64,1280,4096 "
            "with schedule burst"
        ),
        "saturated_nanoseconds_valid": False,
    }


def timing_panel_seed(
    metric: str,
    K: int,
    block_bytes: int,
    schedule: str,
) -> int:
    if (metric not in ("solve", "setup") or
            not isinstance(K, int) or isinstance(K, bool) or
            not K_MIN <= K <= K_MAX or block_bytes not in (64, 1280, 4096) or
            schedule not in SCHEDULES):
        die("timing trace-seed coordinates are outside the frozen grid")
    payload = b"\0".join((
        b"wirehair.wh2.h12-preferred-attempt.timing-trace-seed.v1",
        metric.encode("ascii"), str(K).encode("ascii"),
        str(block_bytes).encode("ascii"), schedule.encode("ascii"),
    ))
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def parse_linux_cpu_list(text: str) -> tuple[int, ...]:
    values: list[int] = []
    for token in text.strip().split(","):
        if not token:
            die("Linux CPU-list contains an empty token")
        fields = token.split("-")
        if len(fields) == 1:
            values.append(strict_uint(fields[0], "Linux CPU-list value"))
        elif len(fields) == 2:
            low = strict_uint(fields[0], "Linux CPU-list low")
            high = strict_uint(fields[1], "Linux CPU-list high")
            if low >= high:
                die("Linux CPU-list range is not increasing")
            values.extend(range(low, high + 1))
        else:
            die("Linux CPU-list token is malformed")
    if tuple(values) != tuple(sorted(set(values))):
        die("Linux CPU-list is duplicate or nonascending")
    return tuple(values)


def format_linux_cpu_list(values: Sequence[int]) -> str:
    cpus = tuple(values)
    if (not cpus or any(
            not isinstance(cpu, int) or isinstance(cpu, bool) or cpu < 0
            for cpu in cpus) or cpus != tuple(sorted(set(cpus)))):
        die("Linux CPU-list values are not canonical")
    ranges: list[str] = []
    first = previous = cpus[0]
    for cpu in cpus[1:]:
        if cpu == previous + 1:
            previous = cpu
            continue
        ranges.append(
            str(first) if first == previous else f"{first}-{previous}")
        first = previous = cpu
    ranges.append(str(first) if first == previous else f"{first}-{previous}")
    return ",".join(ranges)


def parse_systemd_cpu_set(text: str) -> tuple[int, ...]:
    """Parse systemd's canonical space-separated CPUSet serialization."""
    if (not text or text.strip() != text or "," in text or
            "  " in text or any(character.isspace() and character != " "
                                 for character in text)):
        die("systemd CPUSet text is not canonical")
    return parse_linux_cpu_list(text.replace(" ", ","))


def exact_load_filler_pids() -> tuple[int, ...]:
    expected_tail = (b"-c", b"while :; do :; done")
    pids: list[int] = []
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            fields = (entry / "cmdline").read_bytes().rstrip(b"\0").split(b"\0")
            nice = os.getpriority(os.PRIO_PROCESS, int(entry.name))
        except (FileNotFoundError, ProcessLookupError, PermissionError, OSError):
            continue
        if (len(fields) == 3 and Path(os.fsdecode(fields[0])).name == "bash" and
                tuple(fields[1:]) == expected_tail and nice == 19):
            pids.append(int(entry.name))
    return tuple(sorted(pids))


def cpu_busy_ticks(cpu: int) -> int:
    prefix = f"cpu{cpu} "
    try:
        lines = Path("/proc/stat").read_text(encoding="ascii").splitlines()
    except (OSError, UnicodeDecodeError):
        die("cannot read per-CPU scheduler ticks")
    matches = [line for line in lines if line.startswith(prefix)]
    if len(matches) != 1:
        die(f"per-CPU scheduler row is unavailable for CPU {cpu}")
    fields = matches[0].split()[1:]
    if len(fields) < 8 or any(not value.isdigit() for value in fields):
        die("per-CPU scheduler row is malformed")
    counters = [int(value) for value in fields]
    return sum(counters) - counters[3] - counters[4]


class TimingHostSession:
    """Quiesce fillers and reserve one physical core for final timing only."""

    def __init__(
        self,
        contract: dict[str, Any],
        timing: Any,
        core: int,
    ) -> None:
        self.contract = contract
        self.timing = timing
        self.core = core
        self.isolation = timing.inspect_linux_isolation(core)
        self.tools = contract.get("system_tools")
        self.lock_fd: int | None = None
        self.stopped_fillers = 0
        self.filler_cpus: tuple[int, ...] = ()
        self.sibling_states: dict[int, str] = {}
        self.slice_states: dict[str, tuple[str, str]] = {}
        self.previous_signal_handlers: dict[int, Any] = {}
        self.pending_cleanup_signals: list[int] = []
        self.signals_installed = False
        self.closing = False
        self.active = False

    def _handle_timing_signal(self, signum: int, _frame: Any) -> None:
        if self.closing:
            self._defer_timing_signal(signum, _frame)
            return
        if signum == signal.SIGINT:
            raise KeyboardInterrupt
        if signum == signal.SIGTERM:
            die("timing host interrupted by SIGTERM")
        die(f"timing host received unexpected signal {signum}")

    def _defer_timing_signal(self, signum: int, _frame: Any) -> None:
        if signum not in (signal.SIGINT, signal.SIGTERM):
            return
        if signum not in self.pending_cleanup_signals:
            self.pending_cleanup_signals.append(signum)

    def _install_signal_handlers(self) -> None:
        if self.signals_installed or self.previous_signal_handlers:
            die("timing host signal handlers are already installed")
        if threading.current_thread() is not threading.main_thread():
            die("timing host isolation must run on the main thread")
        installed: list[int] = []
        try:
            for signum in (signal.SIGINT, signal.SIGTERM):
                self.previous_signal_handlers[signum] = signal.getsignal(signum)
                signal.signal(signum, self._handle_timing_signal)
                installed.append(signum)
        except BaseException:
            for signum in reversed(installed):
                signal.signal(signum, self.previous_signal_handlers[signum])
            self.previous_signal_handlers.clear()
            raise
        self.signals_installed = True

    def _quiesce_signal_handlers(self) -> None:
        if not self.signals_installed:
            return
        self.closing = True
        errors: list[BaseException] = []
        for signum in self.previous_signal_handlers:
            try:
                # Defer, rather than discard, interruptions while host state
                # is being repaired.  The original handlers are restored and
                # the deferred interruption is surfaced after repair.
                signal.signal(signum, self._defer_timing_signal)
            except BaseException as error:
                errors.append(error)
        if errors:
            raise CampaignError(
                "timing host could not quiesce signal handlers: " +
                "; ".join(str(error) for error in errors))

    def _restore_signal_handlers(self) -> None:
        if not self.signals_installed:
            return
        errors: list[BaseException] = []
        signals = tuple(self.previous_signal_handlers)
        try:
            previous_mask = signal.pthread_sigmask(
                signal.SIG_BLOCK, signals)
        except BaseException as error:
            raise CampaignError(
                f"timing host could not block signals for restoration: {error}"
            ) from error
        failed: dict[int, Any] = {}
        for signum, handler in reversed(tuple(
                self.previous_signal_handlers.items())):
            try:
                signal.signal(signum, handler)
            except BaseException as error:
                errors.append(error)
                failed[signum] = handler
        self.previous_signal_handlers = failed
        self.signals_installed = bool(failed)
        self.closing = bool(failed)
        try:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        except BaseException as error:
            errors.append(error)
        if errors:
            raise CampaignError(
                "timing host could not restore signal handlers: " +
                "; ".join(str(error) for error in errors))

    def tool(self, name: str) -> str:
        if not isinstance(self.tools, dict):
            die("timing host lacks the frozen system-tool inventory")
        record = self.tools.get(name)
        if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
            die(f"timing host lacks frozen {name} identity")
        path = Path(str(record["path"]))
        if (not path.is_absolute() or path.resolve(strict=True) != path or
                common.sha256_file(path) != record["sha256"]):
            die(f"timing host frozen {name} changed")
        return str(path)

    def checked(
        self,
        command: Sequence[str],
        input_bytes: bytes | None = None,
    ) -> bytes:
        result = subprocess.run(
            tuple(command), input=input_bytes, capture_output=True,
            check=False, timeout=60)
        if result.returncode or result.stderr:
            die(
                "timing host command failed: " +
                " ".join(command[:3]) +
                f" (rc={result.returncode})")
        return result.stdout

    def _systemctl_cpu_property(self, unit: str, name: str) -> str:
        if name not in ("AllowedCPUs", "EffectiveCPUs"):
            die("timing requested an unknown systemd CPU property")
        output = self.checked((
            self.tool("systemctl"), "show", f"--property={name}",
            "--value", unit,
        ))
        try:
            value = output.decode("ascii").strip()
        except UnicodeDecodeError:
            die(f"systemd {name} state is not ASCII")
        return value

    def _systemctl_allowed(self, unit: str) -> str:
        return self._systemctl_cpu_property(unit, "AllowedCPUs")

    def _systemctl_effective(self, unit: str) -> str:
        return self._systemctl_cpu_property(unit, "EffectiveCPUs")

    def _set_allowed(self, unit: str, value: str) -> None:
        self.checked((
            self.tool("sudo"), "-n", self.tool("systemctl"),
            "set-property", "--runtime", unit, f"AllowedCPUs={value}",
        ))

    def _write_online(self, cpu: int, value: str) -> None:
        path = f"/sys/devices/system/cpu/cpu{cpu}/online"
        self.checked((
            self.tool("sudo"), "-n", self.tool("tee"), path,
        ), (value + "\n").encode("ascii"))

    def _stop_fillers(self) -> None:
        pids = exact_load_filler_pids()
        expected_cpus = tuple(sorted(os.sched_getaffinity(0)))
        filler_cpus: list[int] = []
        for pid in pids:
            try:
                affinity = tuple(sorted(os.sched_getaffinity(pid)))
            except ProcessLookupError:
                die("timing CPU filler disappeared during topology audit")
            if len(affinity) != 1:
                die("timing CPU filler is not pinned to exactly one CPU")
            filler_cpus.append(affinity[0])
        if pids and tuple(sorted(filler_cpus)) != expected_cpus:
            die("timing CPU fillers do not exactly cover the controller CPU set")
        self.stopped_fillers = len(pids)
        self.filler_cpus = tuple(sorted(filler_cpus))
        self._terminate_pids(pids)
        if exact_load_filler_pids():
            die("timing could not stop the exact CPU filler set")

    def _terminate_pids(self, pids: Sequence[int]) -> None:
        targets = tuple(sorted(set(pids)))
        for pid in targets:
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        deadline = time.monotonic() + 10.0
        while time.monotonic() < deadline:
            for pid in targets:
                try:
                    os.waitpid(pid, os.WNOHANG)
                except (ChildProcessError, ProcessLookupError):
                    pass
            if not any(Path(f"/proc/{pid}").exists() for pid in targets):
                break
            time.sleep(0.02)
        for pid in targets:
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        for pid in targets:
            try:
                os.waitpid(pid, 0)
            except (ChildProcessError, ProcessLookupError):
                pass

    def _restart_fillers(self) -> None:
        if self.stopped_fillers == 0:
            return
        online_cpus = set(parse_linux_cpu_list(
            Path("/sys/devices/system/cpu/online").read_text(
                encoding="ascii").strip()))
        cpus = self.filler_cpus
        if (len(cpus) != self.stopped_fillers or
                not set(cpus).issubset(online_cpus)):
            die("refusing to restart fillers on a changed online CPU set")
        existing = exact_load_filler_pids()
        if existing:
            try:
                existing_cpus = tuple(sorted(
                    next(iter(os.sched_getaffinity(pid)))
                    for pid in existing
                    if len(os.sched_getaffinity(pid)) == 1))
            except ProcessLookupError:
                existing_cpus = ()
            if (len(existing) == self.stopped_fillers and
                    existing_cpus == cpus):
                self.stopped_fillers = 0
                self.filler_cpus = ()
                return
            die("refusing to merge a partial or foreign CPU filler set")
        taskset = self.tool("taskset")
        bash = self.tool("bash")
        last_error = "unknown restart failure"
        for _launch_attempt in range(2):
            launched: list[int] = []
            try:
                for cpu in cpus:
                    process = subprocess.Popen(
                        (taskset, "-c", str(cpu), bash, "-c",
                         "while :; do :; done"),
                        stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL, start_new_session=True)
                    launched.append(process.pid)
                    os.setpriority(os.PRIO_PROCESS, process.pid, 19)
                deadline = time.monotonic() + 5.0
                while (len(exact_load_filler_pids()) != self.stopped_fillers and
                       time.monotonic() < deadline):
                    time.sleep(0.02)
                restarted = exact_load_filler_pids()
                if len(restarted) != self.stopped_fillers:
                    raise CampaignError(
                        "CPU filler set did not restart exactly")
                restarted_cpus = tuple(sorted(
                    next(iter(os.sched_getaffinity(pid))) for pid in restarted
                    if len(os.sched_getaffinity(pid)) == 1))
                if restarted_cpus != cpus:
                    raise CampaignError(
                        "restarted CPU fillers do not cover their original CPUs")
                self.stopped_fillers = 0
                self.filler_cpus = ()
                return
            except (OSError, ProcessLookupError, CampaignError) as error:
                last_error = str(error)
                self._terminate_pids(launched)
                # Reap any launched exact loops that survived long enough to
                # enter the scanner before the direct PID cleanup.
                leftovers = exact_load_filler_pids()
                if leftovers:
                    self._terminate_pids(leftovers)
        die(f"CPU filler restart failed twice: {last_error}")

    def __enter__(self) -> "TimingHostSession":
        if self.active:
            die("timing host session is already active")
        lock_path = Path("/tmp/wirehair-wh2-timing-host.lock")
        self.lock_fd = os.open(
            lock_path, os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0),
            0o600)
        lock_metadata = os.fstat(self.lock_fd)
        if (not stat_module.S_ISREG(lock_metadata.st_mode) or
                lock_metadata.st_nlink != 1):
            os.close(self.lock_fd)
            self.lock_fd = None
            die("timing host lock is not a unique regular file")
        try:
            fcntl.flock(self.lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            os.close(self.lock_fd)
            self.lock_fd = None
            die("another final timing session owns global host isolation")
        try:
            self._install_signal_handlers()
            # This host grants noninteractive sudo per executable rather than
            # through a blanket cached credential, so `sudo -v` is not a
            # meaningful capability probe.  Exercise the exact privileged
            # executables read-only; the mutations below remain fail-closed.
            for name in ("systemctl", "systemd-run", "tee"):
                self.checked((
                    self.tool("sudo"), "-n", self.tool(name), "--version"))
            self._stop_fillers()
            online = set(parse_linux_cpu_list(
                Path("/sys/devices/system/cpu/online").read_text(
                    encoding="ascii").strip()))
            if self.core not in online:
                die("selected timing core is offline")
            # Snapshot cgroup state while every originally-online sibling is
            # still online.  EffectiveCPUs omits offline CPUs, so taking this
            # snapshot later would permanently forget the sibling at restore.
            for unit in ("system.slice", "user.slice", "machine.slice"):
                original = self._systemctl_allowed(unit)
                original_effective = self._systemctl_effective(unit)
                if original_effective:
                    parse_systemd_cpu_set(original_effective)
                self.slice_states[unit] = (original, original_effective)
            for sibling in self.isolation.sibling_cpus:
                path = Path(
                    f"/sys/devices/system/cpu/cpu{sibling}/online")
                state = path.read_text(encoding="ascii").strip()
                if state not in ("0", "1"):
                    die("timing sibling online state is malformed")
                self.sibling_states[sibling] = state
                if state == "1":
                    self._write_online(sibling, "0")
                    if path.read_text(encoding="ascii").strip() != "0":
                        die("timing sibling did not become offline")
            remaining = sorted(online.difference({self.core}).difference(
                self.isolation.sibling_cpus))
            allowed = format_linux_cpu_list(remaining)
            if not allowed:
                die("timing isolation would leave no housekeeping CPU")
            for unit in ("system.slice", "user.slice", "machine.slice"):
                self._set_allowed(unit, allowed)
                observed = self._systemctl_allowed(unit)
                if (not observed or
                        parse_systemd_cpu_set(observed) != tuple(remaining)):
                    die(f"timing failed to exclude its core from {unit}")
                effective = self._systemctl_effective(unit)
                if (unit != "machine.slice" and
                        (not effective or
                         parse_systemd_cpu_set(effective) != tuple(remaining))):
                    die(f"timing effective CPU set changed for {unit}")
            self.active = True
            return self
        except BaseException:
            self.close()
            raise

    def launcher(self) -> tuple[str, ...]:
        if not self.active:
            die("timing host launcher requested outside its session")
        return self.launcher_spec()

    def launcher_spec(self) -> tuple[str, ...]:
        """Return the frozen argv prefix without mutating host isolation."""
        return (
            self.tool("sudo"), "-n", self.tool("systemd-run"),
            "--scope", "--quiet", "--collect",
            "--slice=wirehair-wh2-timing.slice",
            f"--property=AllowedCPUs={self.core}", "--",
            self.tool("numactl"), f"--physcpubind={self.core}",
            f"--membind={self.isolation.numa_node}",
            self.tool("taskset"), "-c", str(self.core),
        )

    def isolation_is_live(self) -> bool:
        if not self.active or exact_load_filler_pids():
            return False
        for sibling in self.isolation.sibling_cpus:
            try:
                if Path(
                        f"/sys/devices/system/cpu/cpu{sibling}/online").read_text(
                            encoding="ascii").strip() != "0":
                    return False
            except OSError:
                return False
        remaining = set(parse_linux_cpu_list(
            Path("/sys/devices/system/cpu/online").read_text(
                encoding="ascii").strip()))
        remaining.discard(self.core)
        remaining.difference_update(self.isolation.sibling_cpus)
        expected = tuple(sorted(remaining))
        for unit in self.slice_states:
            observed = self._systemctl_allowed(unit)
            if not observed or parse_systemd_cpu_set(observed) != expected:
                return False
            effective = self._systemctl_effective(unit)
            if (unit != "machine.slice" and
                    (not effective or
                     parse_systemd_cpu_set(effective) != expected)):
                return False
        return True

    def close(self) -> None:
        errors: list[BaseException] = []
        try:
            self._quiesce_signal_handlers()
        except BaseException as error:
            errors.append(error)
        # Bring siblings online before restoring an EffectiveCPUs snapshot
        # that contains them; otherwise systemd cannot expand the cgroup mask
        # and an apparently blank reset remains permanently narrowed.
        for sibling, state in self.sibling_states.items():
            try:
                current = Path(
                    f"/sys/devices/system/cpu/cpu{sibling}/online").read_text(
                        encoding="ascii").strip()
                if current != state:
                    self._write_online(sibling, state)
                restored = Path(
                    f"/sys/devices/system/cpu/cpu{sibling}/online").read_text(
                        encoding="ascii").strip()
                if restored != state:
                    die("timing sibling online state did not restore")
            except BaseException as error:
                errors.append(error)
        self.sibling_states.clear()
        for unit, (original, original_effective) in reversed(
                tuple(self.slice_states.items())):
            try:
                # systemd clears the AllowedCPUs property text without
                # necessarily expanding cpuset.cpus.effective.  First restore
                # the exact pre-session effective set explicitly; only then
                # restore the original (possibly blank/inherited) property.
                if original_effective:
                    effective_values = parse_systemd_cpu_set(
                        original_effective)
                    self._set_allowed(
                        unit, format_linux_cpu_list(effective_values))
                    if self._systemctl_effective(unit) != original_effective:
                        die(
                            f"timing host did not restore {unit} EffectiveCPUs")
                self._set_allowed(unit, original)
                if self._systemctl_allowed(unit) != original:
                    die(f"timing host did not restore {unit} AllowedCPUs")
                if (original_effective and
                        self._systemctl_effective(unit) != original_effective):
                    die(
                        f"timing host lost restored {unit} EffectiveCPUs")
            except BaseException as error:
                errors.append(error)
        self.slice_states.clear()
        self.active = False
        try:
            self._restart_fillers()
        except BaseException as error:
            errors.append(error)
        if self.lock_fd is not None:
            lock_fd = self.lock_fd
            self.lock_fd = None
            try:
                fcntl.flock(lock_fd, fcntl.LOCK_UN)
            except BaseException as error:
                errors.append(error)
            try:
                os.close(lock_fd)
            except BaseException as error:
                errors.append(error)
        try:
            self._restore_signal_handlers()
        except BaseException as error:
            errors.append(error)
        pending = tuple(self.pending_cleanup_signals)
        self.pending_cleanup_signals.clear()
        if errors:
            interrupted = ""
            if pending:
                interrupted = " (also interrupted by " + ",".join(
                    signal.Signals(signum).name for signum in pending) + ")"
            raise CampaignError(
                "timing host restoration failed: " +
                "; ".join(str(error) for error in errors) + interrupted)
        if pending:
            if pending[0] == signal.SIGINT:
                raise KeyboardInterrupt
            die(
                "timing host interrupted during restoration by " +
                signal.Signals(pending[0]).name)

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> bool:
        self.close()
        return False


@dataclass
class TimingEvidenceToken:
    thermal_mark: dict[str, Any]
    performance_process: subprocess.Popen[bytes]
    sibling_ticks: dict[int, int]
    thermal_path: Path
    performance_path: Path
    performance_command: tuple[str, ...]


class LinuxTimingEvidenceProbe:
    def __init__(
        self,
        timing: Any,
        host: TimingHostSession,
        thermal_log: Path,
        evidence_dir: Path,
    ) -> None:
        self.timing = timing
        self.host = host
        self.thermal_log = thermal_log
        self.evidence_dir = evidence_dir

    def start(
        self,
        spec: Any,
        process_index: int,
        cycle_index: int | None,
    ) -> TimingEvidenceToken:
        if not self.host.isolation_is_live():
            die("timing isolation is not live before an attempt")
        self.evidence_dir.mkdir(parents=True, exist_ok=True)
        suffix = "full" if cycle_index is None else f"cycle-{cycle_index}"
        stem = f"attempt-{process_index:02d}-{suffix}"
        thermal_path = self.evidence_dir / (stem + ".thermal.csv")
        performance_path = self.evidence_dir / (stem + ".performance.bin")
        if thermal_path.exists() or performance_path.exists():
            die("timing evidence path would overwrite an earlier attempt")
        mark = thermal_start(
            self.thermal_log, stale_seconds=5.0, require_zero_edac=False,
            require_zero_dimm_errors=False)
        command = (
            self.host.tool("sudo"), "-n", self.host.tool("turbostat"),
            "--quiet", "--cpu", str(self.host.core), "--interval", "0.1",
            "--show", "Core,CPU,Avg_MHz,Busy%,Bzy_MHz,TSC_MHz",
        )
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True)
        try:
            # Do not prefill the 100 ms frequency window with an idle target
            # core; that would manufacture an APERF/MPERF excursion before the
            # measured process even starts.  A launch-failure check suffices.
            time.sleep(0.02)
            if process.poll() is not None:
                stdout, stderr = process.communicate()
                common.atomic_write(
                    performance_path,
                    b"wirehair-wh2-timing-performance-failed-v1\n" +
                    stdout + b"\n--stderr--\n" + stderr)
                die("turbostat performance monitor exited before timing")
            return TimingEvidenceToken(
                mark, process,
                {cpu: cpu_busy_ticks(cpu)
                 for cpu in self.host.isolation.sibling_cpus},
                thermal_path, performance_path, command)
        except BaseException as error:
            try:
                common.stop_and_reap_process_group(process)
            except BaseException as cleanup_error:
                raise CampaignError(
                    "turbostat cleanup failed after interrupted evidence "
                    f"startup: {cleanup_error}") from error
            raise

    def _performance_samples(
        self,
        raw: bytes,
    ) -> tuple[list[tuple[Decimal, Decimal]], bool]:
        try:
            lines = raw.decode("ascii").splitlines()
        except UnicodeDecodeError:
            return [], True
        header: list[str] | None = None
        samples: list[tuple[Decimal, Decimal]] = []
        gap = False
        required = {"CPU", "Bzy_MHz", "TSC_MHz"}
        for line in lines:
            fields = line.split("\t")
            if required.issubset(fields):
                header = fields
                continue
            if header is None or len(fields) != len(header):
                if line.strip():
                    gap = True
                continue
            row = dict(zip(header, (field.strip() for field in fields)))
            if row.get("CPU") != str(self.host.core):
                continue
            try:
                bzy = Decimal(row["Bzy_MHz"])
                tsc = Decimal(row["TSC_MHz"])
            except (InvalidOperation, KeyError):
                gap = True
                continue
            if not bzy.is_finite() or not tsc.is_finite() or bzy <= 0 or tsc <= 0:
                gap = True
                continue
            samples.append((bzy, tsc))
        if not samples:
            gap = True
        if samples and max(tsc for _bzy, tsc in samples) - min(
                tsc for _bzy, tsc in samples) > Decimal(1):
            gap = True
        return samples, gap

    def finish(self, token: TimingEvidenceToken) -> Any:
        process = token.performance_process
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGINT)
            except ProcessLookupError:
                pass
        try:
            try:
                stdout, stderr = process.communicate(timeout=10)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                stdout, stderr = process.communicate(timeout=5)
        except BaseException as error:
            try:
                common.stop_and_reap_process_group(process)
            except BaseException as cleanup_error:
                raise CampaignError(
                    "turbostat cleanup failed after interrupted evidence "
                    f"collection: {cleanup_error}") from error
            raise
        returncode = process.returncode
        envelope = (
            b"wirehair-wh2-timing-performance-interval-v1\n" +
            f"returncode={returncode}\ncommand_sha256=".encode("ascii") +
            sha256_bytes("\0".join(token.performance_command).encode("utf-8")).encode("ascii") +
            f"\nstdout_bytes={len(stdout)}\nstderr_bytes={len(stderr)}\n".encode("ascii") +
            stdout + b"\n--stderr--\n" + stderr)
        common.atomic_write(token.performance_path, envelope)
        performance_sha256 = common.sha256_file(token.performance_path)
        samples, performance_gap = self._performance_samples(stdout)
        if returncode not in (-signal.SIGINT, 0, 130) or stderr:
            performance_gap = True
        bzy_values = [value for value, _tsc in samples]
        aperf_excursion = bool(
            bzy_values and
            Decimal(100) * max(bzy_values) >
            Decimal(102) * min(bzy_values))
        throttled = any(
            Decimal(100) * bzy < Decimal(95) * tsc
            for bzy, tsc in samples)

        thermal_gap = False
        try:
            # Fast setup/replacement panels can finish between 1 Hz thermal
            # logger samples.  Wait for one post-baseline row so a clean short
            # panel is not deterministically mislabeled as a telemetry gap.
            deadline = time.monotonic() + 2.0
            while (self.thermal_log.stat().st_size <=
                   token.thermal_mark["offset"] and
                   time.monotonic() < deadline):
                time.sleep(0.02)
            summary = thermal_finish(
                self.thermal_log, token.thermal_mark, token.thermal_path,
                limit_c=85.0, stale_seconds=5.0,
                require_zero_edac=False, require_zero_dimm_errors=False,
                dimm_limit_c=90.0)
        except Exception:
            thermal_gap = True
            baseline = token.thermal_mark["baseline"]
            common.atomic_write(
                token.thermal_path,
                token.thermal_mark["header"] +
                token.thermal_mark["baseline_row"])
            summary = {
                "cpu_tctl_max_c": baseline["cpu_tctl_c"],
                "dimm_max_c": max(baseline["dimm_temperatures"]),
                "dimm_read_errors_max": baseline["dimm_read_errors"],
                "edac_ce_delta": 0, "edac_ue_delta": 0,
            }
        thermal_sha256 = common.sha256_file(token.thermal_path)
        edac_ce_delta = int(summary["edac_ce_delta"])
        edac_ue_delta = int(summary["edac_ue_delta"])
        if edac_ce_delta < 0 or edac_ue_delta < 0:
            # Counter rollback/reset is a telemetry discontinuity, never
            # evidence that the interval had no memory-controller errors.
            thermal_gap = True
        sibling_delta = 0
        for cpu, before in token.sibling_ticks.items():
            after = cpu_busy_ticks(cpu)
            if after < before:
                thermal_gap = True
            else:
                sibling_delta += after - before
        live_isolation = self.host.isolation_is_live()
        return self.timing.EnvironmentEvidence(
            core=self.host.core,
            numa_node=self.host.isolation.numa_node,
            exclusive_core=live_isolation,
            load_workers_stopped=not exact_load_filler_pids(),
            sibling_busy_ticks=sibling_delta,
            cpu_temperature_millic=int(round(
                float(summary["cpu_tctl_max_c"]) * 1000.0)),
            dimm_temperature_millic=int(round(
                float(summary["dimm_max_c"]) * 1000.0)),
            dimm_read_errors=int(summary["dimm_read_errors_max"]),
            edac_ce_delta=max(0, edac_ce_delta),
            edac_ue_delta=max(0, edac_ue_delta),
            throttled=throttled,
            aperf_mperf_excursion=aperf_excursion,
            telemetry_gap=thermal_gap or performance_gap,
            thermal_interval_sha256=thermal_sha256,
            performance_interval_sha256=performance_sha256,
        )


def derive_plan(source_root: Path) -> tuple[dict[str, Any], dict[str, bytes]]:
    source_root = source_root.resolve(strict=True)
    verify_source_seals(source_root)
    cohort = derive_active_cohort(source_root)
    active_weights, all_weights = source_timing_weights(source_root, cohort)
    ledgers = serialize_source_ledgers(cohort, active_weights, all_weights)
    protocol = protocol_definition()
    protocol["source"] = {
        "campaign": "sealed-h12-allk-df13c42-v3",
        "commit": SOURCE_COMMIT,
        "failures_sha256": SOURCE_FAILURES_SHA256,
        "paired_cells_sha256": SOURCE_PAIRED_SHA256,
        "phase_seal_sha256": SOURCE_PHASE_SEAL_SHA256,
        "analysis_seal_sha256": SOURCE_ANALYSIS_SEAL_SHA256,
        "contract_sha256": SOURCE_CONTRACT_SHA256,
    }
    protocol["derived_ledgers"] = {
        name: {"bytes": len(data), "sha256": sha256_bytes(data)}
        for name, data in ledgers.items()
    }
    protocol_bytes = common.json_bytes(protocol)
    protocol["self_sha256_excluding_field"] = sha256_bytes(protocol_bytes)
    return protocol, ledgers


def tracked_clean_sources(
    paths: Sequence[Path],
    tool_records: dict[str, dict[str, str]],
) -> tuple[Path, str, str]:
    if not paths:
        die("no tracked source paths")
    git = str(frozen_tool_path(tool_records, "git"))
    repo = Path(subprocess.check_output(
        (git, "-C", str(paths[0].parent), "rev-parse", "--show-toplevel"),
        text=True,
    ).strip()).resolve(strict=True)
    for args, label in (
        (("diff", "--quiet", "HEAD", "--"), "worktree"),
        (("diff", "--cached", "--quiet", "HEAD", "--"), "index"),
    ):
        if subprocess.run((git, "-C", str(repo), *args), check=False).returncode:
            die(f"tracked source {label} is dirty")
    head = subprocess.check_output(
        (git, "-C", str(repo), "rev-parse", "HEAD"), text=True,
    ).strip()
    upstream = subprocess.check_output(
        (git, "-C", str(repo), "rev-parse", "@{upstream}"), text=True,
    ).strip()
    if head != upstream:
        die("prepare requires immutable HEAD already pushed to its upstream")
    branch = subprocess.check_output(
        (git, "-C", str(repo), "symbolic-ref", "--short", "HEAD"),
        text=True,
    ).strip()
    remote = subprocess.check_output(
        (git, "-C", str(repo), "config", "--get", f"branch.{branch}.remote"),
        text=True,
    ).strip()
    merge_ref = subprocess.check_output(
        (git, "-C", str(repo), "config", "--get", f"branch.{branch}.merge"),
        text=True,
    ).strip()
    remote_rows = subprocess.run(
        (git, "-C", str(repo), "ls-remote", "--exit-code", remote, merge_ref),
        text=True, capture_output=True, check=False, timeout=60,
    )
    expected_remote_line = f"{head}\t{merge_ref}\n"
    if remote_rows.returncode or remote_rows.stdout != expected_remote_line:
        die(
            "prepare remote proof failed: upstream branch does not advertise "
            "the exact local HEAD"
        )
    for path in paths:
        resolved = path.resolve(strict=True)
        relative = resolved.relative_to(repo)
        if subprocess.check_output(
                (git, "-C", str(repo), "show", f"HEAD:{relative}")) != \
                common.stable_bytes(resolved):
            die(f"tracked source does not match HEAD: {relative}")
    return repo, head, upstream


def q0_control_identity(binary: Path) -> dict[str, Any]:
    # Compare the dedicated cached-control path against the historical
    # precodefail q0 path.  Timings are ignored; every algebra/work statistic
    # must match exactly at both sides of the adaptive cutoff.
    cases = (4095, 4096, 46252, 64000)
    widths = WIDTHS
    seed = DEV_ROOTS[0]
    schedule = SCHEDULES[0]
    dedicated_header = (
        "N", "bb", "arm", "preferred_attempt", "canonical_attempt",
        "actual_attempt", "routed", "preferred_valid", "fallback", "no_op",
        "direct", "physical_solve", "result", "rank_fail", "error",
        "heavy_shortfall", "inactivated", "binary_def", "heavy_gain",
        "block_xors", "block_muladds",
    )
    legacy_header = (
        "N", "bb", "heavy_family", "mix_count", "overhead", "trials",
        "success", "rank_fail", "error", "fail_rate", "inact_mu",
        "inact_max", "binary_def_mu", "binary_def_max", "heavy_gain_mu",
        "heavy_gain_min", "heavy_shortfall", "solve_ms_mu", "build_ms_mu",
        "peel_ms_mu", "project_ms_mu", "residual_ms_mu", "backsub_ms_mu",
        "seed_attempt", "block_xors_mu", "block_muladds_mu",
        "first_rank_fail", "binary_def_hist", "heavy_gain_hist",
        "failure_trials", "active_packet_peel_seed_xor",
        "mixed_joint_source_xors_mu", "mixed_joint_marginal_xors_mu",
        "mixed_joint_marginal_copies_mu", "mixed_joint_active_deltas_mu",
        "mixed_joint_scratch_bytes_mu", "mixed_dual_source_columns_mu",
    )

    def preamble_tokens(line: str, prefix: str) -> dict[str, str]:
        if not line.startswith(prefix):
            die(f"q0 identity missing {prefix.rstrip()} preamble")
        parsed: dict[str, str] = {}
        for token in line[len(prefix):].split():
            if "=" not in token:
                die(f"q0 identity malformed preamble token {token!r}")
            key, value = token.split("=", 1)
            if not key or not value or key in parsed:
                die(f"q0 identity duplicate/empty preamble token {token!r}")
            parsed[key] = value
        return parsed

    dedicated = subprocess.run([
        str(binary), "preferredattempt", "--mode", "control", "--N",
        ",".join(map(str, cases)), "--bb-list",
        ",".join(map(str, widths)), "--loss", "0.50",
        "--seed", seed, "--schedule", schedule,
    ], capture_output=True, check=False, timeout=600)
    if dedicated.returncode or dedicated.stderr:
        die("dedicated q0 identity command failed")
    dedicated_all_lines = dedicated.stdout.decode("utf-8").splitlines()
    dedicated_preamble = preamble_tokens(
        dedicated_all_lines[0] if dedicated_all_lines else "",
        "# preferredattempt: ")
    dedicated_expected_preamble = {
        "schema": "v2", "policy": "h12-q0-adaptive",
        "canonical": "ascending-first-valid-v1", "mode": "control",
        "loss": "0.5", "seed": seed, "schedule": schedule,
        "preferred_batch_max": "32", "route_cache_sha256": "none",
        "route_context_sha256": "none", "probe_route": "0",
        "physical_solve_accounting": "explicit",
        "systematic_probe_accounting": "explicit",
    }
    if dedicated_preamble != dedicated_expected_preamble:
        die("dedicated q0 identity preamble mismatch")
    dedicated_lines = [
        line for line in dedicated_all_lines
        if not line.startswith("#")
    ]
    dedicated_reader = csv.DictReader(io.StringIO(
        "\n".join(dedicated_lines) + "\n", newline=""))
    if tuple(dedicated_reader.fieldnames or ()) != dedicated_header:
        die("dedicated q0 identity schema mismatch")
    dedicated_rows = list(dedicated_reader)
    if len(dedicated_rows) != len(cases) * len(widths):
        die("dedicated q0 identity cardinality mismatch")
    dedicated_by_key: dict[tuple[int, int], dict[str, str]] = {}
    for row in dedicated_rows:
        key = (
            strict_uint(row["N"], "dedicated identity N"),
            strict_uint(row["bb"], "dedicated identity bb"),
        )
        if key in dedicated_by_key:
            die(f"duplicate dedicated q0 identity row {key}")
        dedicated_by_key[key] = row
    comparisons = 0
    canonical: list[dict[str, Any]] = []
    for K in cases:
        legacy_command = [
            str(binary), "precodefail", "--N", str(K), "--bb-list",
            ",".join(map(str, widths)),
            "--overhead", "0", "--trials", "1", "--threads", "1",
            "--loss", "0.50", "--seed", seed, "--schedule", schedule,
            "--completion", "mixed", "--mix-count", "2",
        ]
        if K >= CUTOFF:
            legacy_command.append("--binary-dense-two-anchor")
        legacy = subprocess.run(
            legacy_command, capture_output=True, check=False, timeout=600)
        if legacy.returncode or legacy.stderr:
            die(f"legacy q0 identity failed at K={K}")
        legacy_all_lines = legacy.stdout.decode("utf-8").splitlines()
        legacy_preamble = preamble_tokens(
            legacy_all_lines[0] if legacy_all_lines else "",
            "# precodefail: ")
        expected_legacy_preamble = {
            "trials": "1", "threads": "1", "loss": "0.5", "seed": seed,
            "completion": "mixed", "mixed_period": "244",
            "mixed_gf256_rows": "10", "mixed_gf16_rows": "2",
            "mixed_geometry": "frozen", "mixed_residue_skew": "0",
            "mixed_residue_schedule": "constant",
            "mixed_residue_hash_seed": "0x0",
            "mixed_residue_hash_keyed": "0",
            "mixed_independent_extension_residues": "0",
            "mixed_residue_buckets_requested": "auto",
            "mixed_extension_residue_seed_xor": "0x4e",
            "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
            "packet_peel_seed_table": "none",
            "binary_dense_rows_override": "0",
            "binary_dense_two_anchor": "1" if K >= CUTOFF else "0",
            "binary_dense_two_anchor_phase": "0",
            "gf256_heavy_rows_override": "0",
            "odd_packet_peel_seed_xor": "0x0",
            "packet_row_seed_multiplier": "0x1",
            "packet_row_seed_avalanche": "0",
            "seed_block_bytes_override": "0", "overhead_stream": "salted",
            "full_payload_solve": "0", "schedule": schedule,
        }
        if legacy_preamble != expected_legacy_preamble:
            die(f"legacy q0 identity preamble mismatch at K={K}")
        lines = [
            line for line in legacy_all_lines
            if not line.startswith("#")
        ]
        legacy_reader = csv.DictReader(io.StringIO(
            "\n".join(lines) + "\n", newline=""))
        if tuple(legacy_reader.fieldnames or ()) != legacy_header:
            die(f"legacy q0 identity schema mismatch at K={K}")
        legacy_rows = list(legacy_reader)
        if len(legacy_rows) != len(widths):
            die(f"legacy q0 identity cardinality mismatch at K={K}")
        legacy_widths: set[int] = set()
        for old in legacy_rows:
            bb = strict_uint(old["bb"], f"legacy identity K={K} bb")
            if bb in legacy_widths:
                die(f"legacy q0 identity duplicate width K={K} bb={bb}")
            legacy_widths.add(bb)
            key = (K, bb)
            row = dedicated_by_key.get(key)
            if row is None or bb not in widths:
                die(f"legacy q0 identity key mismatch {key}")
            fixed_dedicated = {
                "arm": "control", "preferred_attempt": "-1",
                "actual_attempt": row["canonical_attempt"], "routed": "0",
                "preferred_valid": "1", "fallback": "0", "no_op": "0",
                "direct": "0", "physical_solve": "1",
            }
            if any(row.get(field) != expected
                   for field, expected in fixed_dedicated.items()):
                die(f"dedicated q0 control metadata mismatch {key}")
            fixed_legacy = {
                "heavy_family": "periodic", "mix_count": "2",
                "overhead": "0", "trials": "1",
                "active_packet_peel_seed_xor": "0x0",
                "mixed_joint_source_xors_mu": "0.000",
                "mixed_joint_marginal_xors_mu": "0.000",
                "mixed_joint_marginal_copies_mu": "0.000",
                "mixed_joint_active_deltas_mu": "0.000",
                "mixed_joint_scratch_bytes_mu": "0.000",
                "mixed_dual_source_columns_mu": "0.000",
            }
            if any(old.get(field) != expected
                   for field, expected in fixed_legacy.items()):
                die(f"legacy q0 control metadata mismatch {key}")
            mapping = {
                "N": "N", "bb": "bb", "canonical_attempt": "seed_attempt",
                "rank_fail": "rank_fail", "error": "error",
                "heavy_shortfall": "heavy_shortfall",
                "inactivated": "inact_max", "binary_def": "binary_def_max",
                "heavy_gain": "heavy_gain_min", "block_xors": "block_xors_mu",
                "block_muladds": "block_muladds_mu",
            }
            normalized: dict[str, Any] = {"K": K, "BlockBytes": bb}
            for new_field, old_field in mapping.items():
                comparisons += 1
                new_value = row[new_field]
                old_value = old[old_field]
                if old_field in ("block_xors_mu", "block_muladds_mu"):
                    milli = exact_decimal_milli(old_value, old_field)
                    if milli % 1000:
                        die(f"legacy {old_field} is not an exact integer")
                    old_value = str(milli // 1000)
                if new_value != old_value:
                    die(
                        f"q0 control identity mismatch K={K} bb={bb} "
                        f"field={new_field}: {new_value} != {old_value}"
                    )
                normalized[new_field] = new_value
            result = strict_uint(row["result"], "dedicated q0 result")
            rank_fail = strict_uint(row["rank_fail"], "dedicated q0 rank_fail")
            error = strict_uint(row["error"], "dedicated q0 error")
            success = strict_uint(old["success"], "legacy q0 success")
            legacy_rank_fail = strict_uint(
                old["rank_fail"], "legacy q0 rank_fail")
            legacy_error = strict_uint(old["error"], "legacy q0 error")
            if (rank_fail > 1 or error != 0 or success > 1 or
                    legacy_rank_fail > 1 or legacy_error != 0 or
                    success + legacy_rank_fail != 1 or
                    rank_fail != legacy_rank_fail or
                    result != (1 if rank_fail else 0)):
                die(f"q0 identity outcome semantics mismatch {key}")
            inactivated = strict_uint(row["inactivated"], "q0 inactivated")
            binary_def = strict_uint(row["binary_def"], "q0 binary deficit")
            heavy_gain = strict_uint(row["heavy_gain"], "q0 heavy gain")
            aggregate_expected = {
                "fail_rate": "1.00000000" if rank_fail else "0.00000000",
                "inact_mu": f"{inactivated}.000",
                "binary_def_mu": f"{binary_def}.000",
                "heavy_gain_mu": f"{heavy_gain}.000",
                "first_rank_fail": "0" if rank_fail else "-1",
                "binary_def_hist": f"{binary_def}:1",
                "heavy_gain_hist": f"{heavy_gain}:1",
                "failure_trials": "0" if rank_fail else "",
            }
            if any(old.get(field) != expected
                   for field, expected in aggregate_expected.items()):
                die(f"legacy q0 aggregate semantics mismatch {key}")
            canonical.append(normalized)
        if legacy_widths != set(widths):
            die(f"legacy q0 identity width coverage mismatch at K={K}")
    if set(dedicated_by_key) != {
            (K, bb) for K in cases for bb in widths}:
        die("dedicated q0 identity key coverage mismatch")
    return {
        "passed": True, "K": list(cases), "BlockBytes": list(widths),
        "comparisons": comparisons,
        "canonical_sha256": sha256_bytes(common.json_bytes(canonical)),
        "dedicated_stdout_sha256": sha256_bytes(dedicated.stdout),
    }


def command_version(executable: Path, expected: str) -> str:
    result = subprocess.run(
        (str(executable), "--version"), text=True, capture_output=True,
        check=False, timeout=30,
    )
    observed = result.stdout.strip()
    if result.returncode or result.stderr or observed != expected:
        die(
            f"{executable.name} version mismatch: expected {expected!r}, "
            f"observed {observed!r}"
        )
    return observed


def system_tool_identities() -> dict[str, dict[str, str]]:
    tools: dict[str, dict[str, str]] = {}
    for name in (
            "bash", "cmake", "git", "readelf", "ssh", "ssh-add",
            "taskset", "numactl",
            "turbostat", "systemctl", "systemd-run", "sudo", "tee"):
        located = shutil.which(name)
        if located is None:
            die(f"required frozen system tool is unavailable: {name}")
        path = Path(located).resolve(strict=True)
        if not path.is_file() or not os.access(path, os.X_OK):
            die(f"required frozen system tool is not executable: {name}")
        tools[name] = {
            "path": str(path), "sha256": common.sha256_file(path),
        }
    return tools


def frozen_tool_path(
    tool_records: Any,
    name: str,
    *,
    require_current_path: bool = False,
) -> Path:
    if not isinstance(tool_records, dict):
        die("frozen system-tool inventory is missing")
    record = tool_records.get(name)
    if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
        die(f"frozen {name} identity record changed")
    path_text = record.get("path")
    digest = record.get("sha256")
    if (not isinstance(path_text, str) or not isinstance(digest, str) or
            not re.fullmatch(r"[0-9a-f]{64}", digest)):
        die(f"frozen {name} identity is malformed")
    path = Path(path_text)
    if (not path.is_absolute() or path.resolve(strict=True) != path or
            common.sha256_file(path) != digest or not os.access(path, os.X_OK)):
        die(f"frozen {name} executable changed")
    if require_current_path:
        located = shutil.which(name)
        if located is None or Path(located).resolve(strict=True) != path:
            die(f"PATH no longer resolves to the frozen {name} executable")
    return path


def python_runtime_identity() -> dict[str, str]:
    executable = Path(sys.executable).resolve(strict=True)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        die("current Python interpreter is not an executable regular file")
    return {
        "path": str(executable),
        "version": sys.version,
        "sha256": common.sha256_file(executable),
    }


def host_runtime_identity() -> dict[str, Any]:
    boot_id_path = Path("/proc/sys/kernel/random/boot_id")
    try:
        boot_id = boot_id_path.read_text(encoding="ascii").strip()
    except (UnicodeDecodeError, OSError):
        die("host boot ID is not ASCII")
    if not re.fullmatch(
            r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
            boot_id):
        die("host boot ID is malformed")
    uname = os.uname()
    return {
        "boot_id": boot_id,
        "uname": {
            "sysname": uname.sysname, "nodename": uname.nodename,
            "release": uname.release, "version": uname.version,
            "machine": uname.machine,
        },
    }


def obtain_drand_verifier(
    wrapper: Path,
    package_json: Path,
    package_lock: Path,
) -> dict[str, Any]:
    """Install the exact lock, hash its CJS bundle, and run the offline KAT."""
    node_name = shutil.which("node")
    npm_name = shutil.which("npm")
    if node_name is None or npm_name is None:
        die("pinned drand verification requires node and npm")
    node = Path(node_name).resolve(strict=True)
    npm = Path(npm_name).resolve(strict=True)
    node_version = command_version(node, DRAND_NODE_VERSION)
    npm_version = command_version(npm, DRAND_NPM_VERSION)
    with tempfile.TemporaryDirectory(prefix="wh2-drand-client-") as temporary:
        temp = Path(temporary)
        common.atomic_write(temp / "package.json", common.stable_bytes(package_json))
        common.atomic_write(
            temp / "package-lock.json", common.stable_bytes(package_lock))
        result = subprocess.run(
            (
                str(npm), "ci", "--ignore-scripts", "--no-audit", "--no-fund",
            ),
            cwd=temp, text=True, capture_output=True, check=False, timeout=120,
        )
        if result.returncode:
            die(f"failed to install pinned drand-client lock: {result.stderr.strip()}")
        installed_package = temp / "node_modules/drand-client/package.json"
        try:
            installed = json.loads(common.stable_bytes(installed_package))
        except json.JSONDecodeError:
            die("installed drand-client package metadata is malformed")
        if (not isinstance(installed, dict) or
                installed.get("name") != "drand-client" or
                installed.get("version") != DRAND_CLIENT_VERSION):
            die("installed drand-client package version mismatch")
        bundle_path = (temp / "node_modules/drand-client/build/cjs/index.cjs")
        if bundle_path.is_symlink():
            die("installed drand-client CJS bundle is a symlink")
        bundle_path = bundle_path.resolve(strict=True)
        if temp.resolve(strict=True) not in bundle_path.parents:
            die("installed drand-client CJS bundle escaped its install root")
        bundle = bundle_path.read_bytes()
        if (len(bundle) != DRAND_CLIENT_BUNDLE_BYTES or
                sha256_bytes(bundle) != DRAND_CLIENT_BUNDLE_SHA256):
            die("installed drand-client CJS bundle content mismatch")
        selftest = subprocess.run(
            (str(node), str(wrapper), str(bundle_path), "offline-selftest"),
            text=True, capture_output=True, check=False, timeout=90,
        )
        if selftest.returncode or selftest.stderr:
            die(f"pinned drand offline self-test failed: {selftest.stderr.strip()}")
        try:
            selftest_record = json.loads(selftest.stdout)
        except json.JSONDecodeError:
            die("pinned drand offline self-test returned malformed JSON")
        if (not isinstance(selftest_record, dict) or
                selftest_record.get("schema") !=
                "wirehair.wh2.drand_quicknet_offline_kat.v1" or
                selftest_record.get("beacon") != DRAND_KNOWN_ROUND or
                selftest_record.get("canonical_sha256") !=
                DRAND_CANONICAL_KNOWN_SHA256):
            die("pinned drand offline self-test record mismatch")
    return {
        "bundle": bundle,
        "node": str(node),
        "node_version": node_version,
        "node_sha256": common.sha256_file(node),
        "npm": str(npm),
        "npm_version": npm_version,
        "npm_sha256": common.sha256_file(npm),
        "offline_selftest": selftest_record,
    }


def fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    descriptor = os.open(path, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_tree(root: Path) -> None:
    """Flush every staged file and directory before publishing the tree."""
    directories = [root]
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            die(f"refusing symlink in staged result: {path}")
        if path.is_dir():
            directories.append(path)
            continue
        if not path.is_file():
            die(f"refusing non-file in staged result: {path}")
        descriptor = os.open(path, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    for directory in sorted(
            directories, key=lambda item: len(item.parts), reverse=True):
        fsync_directory(directory)


def fixed_json_bytes(value: Any) -> bytes:
    """Canonical compact JSON whose insertion order is part of the schema."""
    return (
        json.dumps(
            value, separators=(",", ":"), allow_nan=False,
            ensure_ascii=True,
        ) + "\n"
    ).encode("ascii")


def load_json_object(path: Path, context: str) -> dict[str, Any]:
    try:
        value = json.loads(common.stable_bytes(path))
    except (json.JSONDecodeError, UnicodeDecodeError, OSError):
        die(f"{context} is not valid JSON")
    if not isinstance(value, dict):
        die(f"{context} must be a JSON object")
    return value


@contextmanager
def holdout_controller_lock(
    result_dir: Path,
    phase: str,
) -> Iterator[None]:
    beacon_dir = result_dir / "beacon"
    beacon_dir.mkdir(mode=0o700, exist_ok=True)
    if beacon_dir.is_symlink():
        die("beacon controller directory must not be a symlink")
    lock_path = beacon_dir / f"{phase}.lock"
    descriptor = os.open(
        lock_path, os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0), 0o600)
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            die(f"another {phase} holdout controller is active")
        metadata = os.fstat(descriptor)
        if (not stat_module.S_ISREG(metadata.st_mode) or
                metadata.st_nlink != 1 or metadata.st_size != 0):
            die("holdout controller lock is not a unique empty regular file")
        yield
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def is_relative_prefix(relative: str, prefix: str) -> bool:
    return relative == prefix or relative.startswith(prefix + "/")


def validate_manifest_path(relative: str) -> bytes:
    if (not relative or relative.startswith("/") or "\\" in relative or
            any(character in relative for character in "\r\n\0") or
            unicodedata.normalize("NFC", relative) != relative):
        die(f"noncanonical phase-manifest path {relative!r}")
    parts = relative.split("/")
    if any(part in ("", ".", "..") for part in parts):
        die(f"unsafe phase-manifest path {relative!r}")
    return relative.encode("utf-8")


def phase_manifest_bytes(
    result_dir: Path,
    spec: HoldoutPhaseSpec,
) -> bytes:
    allowed_top = {
        "prepare.json", "frozen", "development", "jobs", "ledgers",
        "completed", "thermal", "beacon", "seals", "holdout",
    }
    if spec.phase == "h2":
        allowed_top.add("h1")
    records: list[tuple[bytes, str, int, str]] = []
    seen: set[str] = set()
    for path in sorted(result_dir.rglob("*")):
        relative = path.relative_to(result_dir).as_posix()
        relative_bytes = validate_manifest_path(relative)
        before = path.lstat()
        if path.is_symlink():
            die(f"symlink is forbidden in phase inputs: {relative}")
        if any(is_relative_prefix(relative, prefix)
               for prefix in spec.excluded_prefixes):
            continue
        top = relative.split("/", 1)[0]
        if top not in allowed_top:
            die(f"unexpected top-level phase input {relative!r}")
        if path.is_dir():
            continue
        if not path.is_file():
            die(f"nonregular phase input is forbidden: {relative}")
        data = common.stable_bytes(path)
        after = path.lstat()
        identity_before = (
            before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
        )
        identity_after = (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
        )
        if identity_before != identity_after or len(data) != after.st_size:
            die(f"phase input changed while hashing: {relative}")
        if relative in seen:
            die(f"duplicate phase input path: {relative}")
        seen.add(relative)
        records.append((
            relative_bytes, sha256_bytes(data), len(data), relative,
        ))
    missing = sorted(set(spec.required_files).difference(seen))
    if missing:
        die(f"phase inputs are missing required files: {', '.join(missing)}")
    records.sort(key=lambda row: row[0])
    return MANIFEST_MAGIC + b"".join(
        f"{digest} {size} ".encode("ascii") + relative_bytes + b"\n"
        for relative_bytes, digest, size, _relative in records
    )


def verify_manifest_bytes(
    result_dir: Path,
    spec: HoldoutPhaseSpec,
    expected: bytes,
) -> None:
    if phase_manifest_bytes(result_dir, spec) != expected:
        die(f"PERMANENT: {spec.phase} sealed input manifest changed")


def verify_frozen_staged_anchor(
    result_dir: Path,
    prepare_record: dict[str, Any],
) -> None:
    frozen = result_dir / "frozen"
    if frozen.is_symlink() or not frozen.is_dir():
        die("PERMANENT: frozen staged directory identity changed")
    staged_path = frozen / "staged.sha256"
    staged_entries = common.verify_sha_manifest(result_dir, staged_path)
    if (common.sha256_file(staged_path) !=
            prepare_record.get("staged_seal_sha256")):
        die("PERMANENT: frozen staged-seal trust anchor changed")
    staged_paths = [
        path for path in frozen.iterdir()
        if path.name not in ("staged.sha256", "r1_freeze_publication.json")
    ]
    if any(path.is_symlink() or not path.is_file() for path in staged_paths):
        die("PERMANENT: frozen staged file identity changed")
    expected_staged = {path.resolve(strict=True) for path in staged_paths}
    if set(staged_entries) != expected_staged:
        die("PERMANENT: frozen staged-seal file coverage changed")
    if (common.sha256_file(frozen / "contract.json") !=
            prepare_record.get("contract_sha256")):
        die("PERMANENT: frozen contract trust anchor changed")


def verify_prepare_contract_bindings(
    frozen: Path,
    prepare_record: dict[str, Any],
    contract: dict[str, Any],
) -> None:
    duplicated_bindings = {
        "source_commit": contract.get("source_commit"),
        "binary_sha256": contract.get("binary_sha256"),
        "script_sha256": contract.get("script_sha256"),
        "allk_sha256": contract.get("allk_sha256"),
        "helper_sha256": contract.get("helper_sha256"),
        "campaign_module_sha256": contract.get("campaign_module_sha256"),
        "holdout_module_sha256": contract.get("holdout_module_sha256"),
        "timing_module_sha256": contract.get("timing_module_sha256"),
        "protocol_sha256": contract.get("protocol_sha256"),
        "route_context_sha256": contract.get("route_context_sha256"),
        "q0_identity_sha256": contract.get("q0_identity_sha256"),
        "drand_wrapper_sha256": contract.get("drand_wrapper_sha256"),
        "drand_package_sha256": contract.get("drand_package_sha256"),
        "drand_package_lock_sha256":
            contract.get("drand_package_lock_sha256"),
        "drand_client_bundle_sha256":
            contract.get("drand_client_bundle_sha256"),
        "node": contract.get("node"),
        "npm": contract.get("npm"),
        "python": contract.get("python"),
        "system_tools": contract.get("system_tools"),
        "host_runtime": contract.get("host_runtime"),
        "drand_offline_selftest": contract.get("drand_offline_selftest"),
    }
    if any(prepare_record.get(key) != value
           for key, value in duplicated_bindings.items()):
        die("PERMANENT: prepare/contract identity binding changed")
    q0_identity = load_json_object(frozen / "q0_identity.json", "q0 identity")
    if (q0_identity != prepare_record.get("q0_identity") or
            common.sha256_file(frozen / "q0_identity.json") !=
                prepare_record.get("q0_identity_sha256")):
        die("PERMANENT: q0 identity trust anchor changed")


def verify_frozen_controller_runtime(
    result_dir: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    frozen = result_dir / "frozen"
    prepare_record = load_json_object(result_dir / "prepare.json", "prepare record")
    contract = load_json_object(frozen / "contract.json", "frozen contract")
    verify_frozen_staged_anchor(result_dir, prepare_record)
    exact_hashes = {
        frozen / "wirehair_v2_bench": contract.get("binary_sha256"),
        frozen / "wh2_preferred_attempt_search.py": contract.get("script_sha256"),
        frozen / "wh2_preferred_attempt_campaign.py":
            contract.get("campaign_module_sha256"),
        frozen / "wh2_preferred_attempt_holdout.py":
            contract.get("holdout_module_sha256"),
        frozen / "wh2_preferred_attempt_timing.py":
            contract.get("timing_module_sha256"),
        frozen / "wh2_drand_verify.cjs": contract.get("drand_wrapper_sha256"),
        frozen / "drand-client.cjs": contract.get("drand_client_bundle_sha256"),
        frozen / "drand-package.json": contract.get("drand_package_sha256"),
        frozen / "drand-package-lock.json":
            contract.get("drand_package_lock_sha256"),
    }
    for path, expected in exact_hashes.items():
        if (not isinstance(expected, str) or
                common.sha256_file(path.resolve(strict=True)) != expected):
            die(f"PERMANENT: frozen controller input changed: {path.name}")
    verify_prepare_contract_bindings(frozen, prepare_record, contract)
    current_script = Path(__file__).resolve(strict=True)
    expected_script = (
        frozen / "wh2_preferred_attempt_search.py").resolve(strict=True)
    if (current_script != expected_script or
            common.sha256_file(current_script) != contract.get("script_sha256")):
        die("PERMANENT: controller must execute the exact frozen script")
    node_record = contract.get("node")
    npm_record = contract.get("npm")
    if not isinstance(node_record, dict) or not isinstance(npm_record, dict):
        die("frozen runtime records are missing")
    for record, expected_version in (
            (node_record, DRAND_NODE_VERSION), (npm_record, DRAND_NPM_VERSION)):
        executable = Path(str(record.get("path", ""))).resolve(strict=True)
        if (record.get("version") != expected_version or
                common.sha256_file(executable) != record.get("sha256") or
                command_version(executable, expected_version) != expected_version):
            die(f"PERMANENT: frozen {executable.name} runtime changed")
    tool_records = contract.get("system_tools")
    expected_tool_names = {
        "bash", "cmake", "git", "readelf", "ssh", "ssh-add",
        "taskset", "numactl",
        "turbostat", "systemctl", "systemd-run", "sudo", "tee",
    }
    if not isinstance(tool_records, dict) or set(tool_records) != expected_tool_names:
        die("PERMANENT: frozen system-tool inventory changed")
    for name in sorted(expected_tool_names):
        frozen_tool_path(tool_records, name)
    python_record = contract.get("python")
    if (not isinstance(python_record, dict) or
            set(python_record) != {"path", "version", "sha256"} or
            python_record != python_runtime_identity()):
        die("PERMANENT: frozen Python interpreter identity changed")
    if contract.get("host_runtime") != host_runtime_identity():
        die("PERMANENT: frozen host boot/kernel identity changed")
    if (prepare_record.get("holdout_root_files_present") is not False or
            contract.get("holdout_root_files_present") is not False):
        die("prepare/frozen contract binding mismatch")
    verify_r1_freeze_publication(result_dir, prepare_record, contract)
    return prepare_record, contract


def load_frozen_python_module(
    result_dir: Path,
    module_name: str,
    contract_key: str,
) -> Any:
    """Import one staged helper only from the verified frozen directory."""
    _prepare, contract = verify_frozen_controller_runtime(result_dir)
    if not re.fullmatch(r"wh2_preferred_attempt_[a-z_]+", module_name):
        die("frozen helper module name is not canonical")
    expected = (result_dir / "frozen" / (module_name + ".py")).resolve(
        strict=True)
    expected_hash = contract.get(contract_key)
    if (not isinstance(expected_hash, str) or
            not re.fullmatch(r"[0-9a-f]{64}", expected_hash) or
            common.sha256_file(expected) != expected_hash):
        die(f"PERMANENT: frozen helper binding changed: {module_name}")
    module = importlib.import_module(module_name)
    module_file = getattr(module, "__file__", None)
    if (not isinstance(module_file, str) or
            Path(module_file).resolve(strict=True) != expected or
            common.sha256_file(expected) != expected_hash):
        die(f"PERMANENT: imported helper is not frozen: {module_name}")
    return module


def load_frozen_cohort_bins(
    result_dir: Path,
    campaign: Any,
) -> tuple[tuple[int, ...], dict[int, tuple[int, ...]]]:
    frozen = result_dir / "frozen"
    cohort_bytes = common.stable_bytes(frozen / "cohort.txt")
    bins_bytes = common.stable_bytes(frozen / "bins.tsv")
    if (sha256_bytes(cohort_bytes) != COHORT_SHA256 or
            sha256_bytes(bins_bytes) != BINS_SHA256):
        die("PERMANENT: frozen development cohort/bin ledger changed")
    try:
        cohort_lines = cohort_bytes.decode("ascii").splitlines()
        cohort = tuple(
            strict_uint(value, "frozen cohort K") for value in cohort_lines)
    except UnicodeDecodeError:
        die("frozen cohort ledger is not ASCII")
    if (len(cohort) != ACTIVE_K_COUNT or tuple(sorted(cohort)) != cohort or
            len(set(cohort)) != len(cohort) or
            any(not CUTOFF <= K <= K_MAX for K in cohort)):
        die("frozen development cohort geometry changed")
    bins: dict[int, tuple[int, ...]] = {}
    try:
        bin_lines = bins_bytes.decode("ascii").splitlines()
    except UnicodeDecodeError:
        die("frozen bin ledger is not ASCII")
    if len(bin_lines) != BIN_COUNT:
        die("frozen bin ledger cardinality changed")
    for number, line in enumerate(bin_lines):
        fields = line.split("\t")
        if len(fields) != 4:
            die("frozen bin ledger row schema changed")
        index = strict_uint(fields[0], "frozen bin index")
        _load = strict_uint(fields[1], "frozen bin load")
        count = strict_uint(fields[2], "frozen bin count")
        values = tuple(
            strict_uint(value, "frozen bin K")
            for value in fields[3].split(",")
        )
        if index != number or count != len(values):
            die("frozen bin ledger index/count changed")
        bins[index] = values
    campaign.validate_bins(bins, cohort)
    return cohort, bins


def load_frozen_allk_bins(
    result_dir: Path,
    campaign: Any,
) -> tuple[tuple[int, ...], dict[int, tuple[int, ...]]]:
    frozen = result_dir / "frozen"
    data = common.stable_bytes(frozen / "allk_bins.tsv")
    if sha256_bytes(data) != ALLK_BINS_SHA256:
        die("PERMANENT: frozen all-K bin ledger changed")
    try:
        lines = data.decode("ascii").splitlines()
    except UnicodeDecodeError:
        die("frozen all-K bin ledger is not ASCII")
    if len(lines) != BIN_COUNT:
        die("frozen all-K bin ledger cardinality changed")
    bins: dict[int, tuple[int, ...]] = {}
    for number, line in enumerate(lines):
        fields = line.split("\t")
        if len(fields) != 4:
            die("frozen all-K bin ledger row schema changed")
        index = strict_uint(fields[0], "frozen all-K bin index")
        _load = strict_uint(fields[1], "frozen all-K bin load")
        count = strict_uint(fields[2], "frozen all-K bin count")
        values = tuple(
            strict_uint(value, "frozen all-K bin K")
            for value in fields[3].split(","))
        if index != number or count != len(values):
            die("frozen all-K bin ledger index/count changed")
        bins[index] = values
    domain = tuple(range(K_MIN, K_MAX + 1))
    campaign.validate_bins(bins, domain)
    return domain, bins


def now_utc_ms() -> int:
    return time.time_ns() // 1_000_000


def clean_source_commit(contract: dict[str, Any]) -> tuple[Path, str]:
    repo = Path(str(contract.get("source_repo", ""))).resolve(strict=True)
    commit = str(contract.get("source_commit", ""))
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        die("frozen source commit is invalid")
    git = str(frozen_tool_path(contract.get("system_tools"), "git"))
    head = subprocess.check_output(
        (git, "-C", str(repo), "rev-parse", "HEAD"), text=True,
    ).strip()
    if head != commit:
        die("PERMANENT: source worktree no longer names the frozen commit")
    for arguments in (("diff", "--quiet", "HEAD", "--"),
                      ("diff", "--cached", "--quiet", "HEAD", "--")):
        if subprocess.run(
                (git, "-C", str(repo), *arguments), check=False).returncode:
            die("PERMANENT: source worktree or index became dirty")
    return repo, commit


def github_remote_url(repo: Path, git_path: Path | str) -> str:
    url = subprocess.check_output(
        (str(git_path), "-C", str(repo), "remote", "get-url", "origin"),
        text=True,
    ).strip()
    if not (re.fullmatch(r"git@github\.com:[^\s]+", url) or
            re.fullmatch(r"ssh://git@github\.com/[^\s]+", url)):
        die("holdout seal publication requires origin over github.com SSH")
    return url


def build_seal_record(
    phase: str,
    manifest_sha256: str,
    contract: dict[str, Any],
    remote_url: str,
    seal_ms: int,
    latest_bound: dict[str, Any],
    abandoned_round: int = 0,
) -> tuple[dict[str, Any], bytes]:
    not_before_ms = seal_ms + 120_000
    local_formula_round = target_beacon_round(seal_ms)
    latest_round = latest_bound.get("latest_round")
    if (not isinstance(latest_round, int) or isinstance(latest_round, bool) or
            not 1 <= latest_round <= (1 << 53) - 41 or
            not isinstance(abandoned_round, int) or
            isinstance(abandoned_round, bool) or abandoned_round < 0):
        die("verified-latest or abandoned round is invalid")
    latest_minimum_round = latest_round + 40
    round_number = max(
        local_formula_round, latest_minimum_round, abandoned_round + 1)
    round_time_ms = (
        DRAND_GENESIS_TIME * 1000 +
        (round_number - 1) * DRAND_PERIOD_SECONDS * 1000
    )
    record = {
        "schema": "wirehair.wh2.holdout_seal.v1",
        "phase": phase,
        "manifest_sha256": manifest_sha256,
        "source_commit": contract["source_commit"],
        "dirty_diff_sha256": EMPTY_SHA256,
        "protocol_sha256": contract["protocol_sha256"],
        "binary_sha256": contract["binary_sha256"],
        "drand_wrapper_sha256": contract["drand_wrapper_sha256"],
        "drand_client_bundle_sha256": contract["drand_client_bundle_sha256"],
        "drand_package_sha256": contract["drand_package_sha256"],
        "drand_package_lock_sha256": contract["drand_package_lock_sha256"],
        "node_version": contract["node"]["version"],
        "npm_version": contract["npm"]["version"],
        "drand_client_sri": DRAND_CLIENT_INTEGRITY,
        "chain_info": {
            "public_key": DRAND_PUBLIC_KEY,
            "period": DRAND_PERIOD_SECONDS,
            "genesis_time": DRAND_GENESIS_TIME,
            "hash": DRAND_CHAIN_HASH,
            "groupHash": DRAND_GROUP_HASH,
            "schemeID": "bls-unchained-g1-rfc9380",
            "metadata": {"beaconID": "quicknet"},
        },
        "origins": [f"https://{host}" for host in DRAND_RELAYS],
        "github_remote": remote_url,
        "seal_ms": seal_ms,
        "not_before_ms": not_before_ms,
        "local_formula_round": local_formula_round,
        "verified_latest": latest_bound,
        "verified_latest_minimum_round": latest_minimum_round,
        "abandoned_round_exclusive": abandoned_round,
        "round": round_number,
        "round_time_ms": round_time_ms,
    }
    encoded = fixed_json_bytes(record)
    decoded = json.loads(encoded)
    if decoded != record or fixed_json_bytes(decoded) != encoded:
        die("internal fixed-order holdout seal encoding mismatch")
    if (round_number < local_formula_round or
            round_number < latest_minimum_round or
            round_number <= latest_round or
            round_number > (1 << 53) - 1 or
            (abandoned_round and round_number <= abandoned_round) or
            round_time_ms < not_before_ms):
        die("internal holdout seal target-time mismatch")
    return record, encoded


def obtain_verified_latest_bound(
    contract: dict[str, Any],
    result_dir: Path,
) -> dict[str, Any]:
    """Verify latest only as a lower bound; never return its beacon entropy."""
    command = (
        str(contract["node"]["path"]),
        str(result_dir / "frozen/wh2_drand_verify.cjs"),
        str(result_dir / "frozen/drand-client.cjs"),
        "latest-bound",
    )
    try:
        result = subprocess.run(
            command, capture_output=True, check=False, timeout=30)
    except subprocess.TimeoutExpired:
        die("verified-latest Quicknet quorum timed out; no seal was created")
    try:
        wave = json.loads(result.stdout)
    except (json.JSONDecodeError, UnicodeDecodeError):
        die("verified-latest Quicknet wrapper returned malformed JSON")
    if result.returncode or result.stderr or not isinstance(wave, dict):
        die("verified-latest Quicknet quorum unavailable; no seal was created")
    beacon = wave.get("beacon")
    consensus = wave.get("consensus_origins")
    observations = wave.get("observations")
    all_origins = [f"https://{host}" for host in DRAND_RELAYS]
    expected_wave_keys = {
        "schema", "role", "chain_hash", "status", "latest_round",
        "consensus_origins", "observations", "beacon",
    }
    if (set(wave) != expected_wave_keys or
            wave.get("status") != "QUORUM" or
            wave.get("schema") !=
                "wirehair.wh2.drand_quicknet_latest_wave.v1" or
            wave.get("role") !=
                "preseal-clock-liveness-lower-bound-only" or
            wave.get("chain_hash") != DRAND_CHAIN_HASH or
            not isinstance(beacon, dict) or
            wave.get("latest_round") != beacon.get("round") or
            not isinstance(consensus, list) or len(consensus) < 3 or
            len(consensus) != len(set(consensus)) or
            consensus != [origin for origin in all_origins
                           if origin in consensus] or
            not isinstance(observations, list) or
            len(observations) != len(all_origins) or
            {item.get("origin") for item in observations
             if isinstance(item, dict)} != set(all_origins)):
        die("verified-latest Quicknet result has invalid quorum coverage")
    canonical = canonical_beacon_bytes(beacon)
    offline_verify_beacon(contract, result_dir, beacon)
    raw_hashes: dict[str, str] = {}
    by_origin = {
        item["origin"]: item for item in observations if isinstance(item, dict)
    }
    for origin in all_origins:
        item = by_origin[origin]
        digest = item.get("raw_response_sha256")
        if origin in consensus and (
                item.get("ok") is not True or item.get("beacon") != beacon or
                item.get("canonical_sha256") != sha256_bytes(canonical)):
            die("verified-latest quorum observation disagrees with consensus")
        if origin in consensus and (
                not isinstance(digest, str) or
                not re.fullmatch(r"[0-9a-f]{64}", digest)):
            die("verified-latest quorum lacks a raw response hash")
        if isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest):
            raw_hashes[origin] = digest
    return {
        "schema": "wirehair.wh2.verified_latest_bound.v1",
        "latest_round": beacon["round"],
        "canonical_beacon_sha256": sha256_bytes(canonical),
        "consensus_origins": [
            origin for origin in all_origins if origin in set(consensus)
        ],
        "raw_response_sha256": raw_hashes,
        "wrapper_result_sha256": sha256_bytes(result.stdout),
        "verified_ms": now_utc_ms(),
    }


def seal_tag_annotation(
    phase: str,
    seal_sha256: str,
    manifest_sha256: str,
    seal_record: dict[str, Any],
) -> str:
    return (
        "wirehair-wh2-holdout-seal-tag-v1\n"
        f"phase={phase}\n"
        f"seal_record_sha256={seal_sha256}\n"
        f"manifest_sha256={manifest_sha256}\n"
        f"round={seal_record['round']}\n"
        f"round_time_ms={seal_record['round_time_ms']}\n"
        f"source_commit={seal_record['source_commit']}\n"
    )


def newest_github_agent_environment(
    tool_records: Any,
) -> dict[str, str]:
    # A child shell cannot update this process's environment.  Ignore the
    # usually-stale persisted SSH_AUTH_SOCK and scan live forwarded sockets.
    ssh = str(frozen_tool_path(tool_records, "ssh"))
    ssh_add = str(frozen_tool_path(tool_records, "ssh-add"))
    candidates: list[tuple[int, Path]] = []
    for path in Path("/tmp").glob("ssh-*/agent.*"):
        try:
            if path.is_socket():
                candidates.append((path.stat().st_mtime_ns, path))
        except OSError:
            continue
    candidates.sort(key=lambda item: item[0], reverse=True)
    for _mtime_ns, socket_path in candidates:
        environment = os.environ.copy()
        environment["SSH_AUTH_SOCK"] = str(socket_path)
        # Override user/repository Git SSH configuration with the exact frozen
        # client.  Its absolute path does not depend on the ambient PATH.
        environment.pop("GIT_SSH", None)
        environment["GIT_SSH_COMMAND"] = shlex.join((
            ssh, "-o", "BatchMode=yes", "-o", "ConnectTimeout=10"))
        environment["GIT_SSH_VARIANT"] = "ssh"
        listed = subprocess.run(
            (ssh_add, "-l"), env=environment, text=True,
            capture_output=True, check=False, timeout=15,
        )
        if listed.returncode:
            continue
        authenticated = subprocess.run(
            (
                ssh, "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                "-T", "git@github.com",
            ),
            env=environment, text=True, capture_output=True,
            check=False, timeout=20,
        )
        combined = authenticated.stdout + authenticated.stderr
        if (authenticated.returncode == 1 and
                "successfully authenticated" in combined):
            return environment
    die("no current forwarded SSH agent authenticates to GitHub")


def tag_remote_rows(
    repo: Path,
    tag_name: str,
    environment: dict[str, str],
    tool_records: Any,
) -> tuple[bytes, dict[str, str]]:
    git = str(frozen_tool_path(tool_records, "git"))
    result = subprocess.run(
        (
            git, "-C", str(repo), "ls-remote", "--tags", "origin",
            f"refs/tags/{tag_name}", f"refs/tags/{tag_name}^{{}}",
        ),
        env=environment, capture_output=True, check=False, timeout=60,
    )
    if result.returncode or result.stderr:
        die("GitHub tag verification failed")
    rows: dict[str, str] = {}
    for line in result.stdout.decode("ascii").splitlines():
        fields = line.split("\t")
        if (len(fields) != 2 or not re.fullmatch(r"[0-9a-f]{40}", fields[0]) or
                fields[1] in rows):
            die("GitHub tag verification returned noncanonical rows")
        rows[fields[1]] = fields[0]
    return result.stdout, rows


def publish_seal_tag(
    repo: Path,
    phase: str,
    seal_sha256: str,
    manifest_sha256: str,
    seal_record: dict[str, Any],
    tool_records: Any,
) -> dict[str, Any]:
    git = str(frozen_tool_path(tool_records, "git"))
    environment = newest_github_agent_environment(tool_records)
    tag_name = f"wh2-h12-kpreferred-{phase}-seal-{seal_sha256[:16]}"
    if not re.fullmatch(r"wh2-h12-kpreferred-(?:h1|h2)-seal-[0-9a-f]{16}", tag_name):
        die("internal holdout seal tag name mismatch")
    annotation = seal_tag_annotation(
        phase, seal_sha256, manifest_sha256, seal_record)
    local_ref = f"refs/tags/{tag_name}"
    local_exists = subprocess.run(
        (git, "-C", str(repo), "show-ref", "--verify", "--quiet", local_ref),
        check=False,
    ).returncode == 0
    if not local_exists:
        subprocess.run(
            (
                git, "-C", str(repo), "tag", "-a", tag_name,
                seal_record["source_commit"], "-m", annotation,
            ),
            check=True, timeout=30,
        )
    tag_object = subprocess.check_output(
        (git, "-C", str(repo), "rev-parse", f"{tag_name}^{{tag}}"),
        text=True,
    ).strip()
    peeled = subprocess.check_output(
        (git, "-C", str(repo), "rev-parse", f"{tag_name}^{{}}"),
        text=True,
    ).strip()
    tag_bytes = subprocess.check_output(
        (git, "-C", str(repo), "cat-file", "tag", tag_object))
    try:
        tag_message = tag_bytes.split(b"\n\n", 1)[1].decode("utf-8")
    except (IndexError, UnicodeDecodeError):
        die("local annotated seal tag is malformed")
    if (peeled != seal_record["source_commit"] or tag_message != annotation or
            not re.fullmatch(r"[0-9a-f]{40}", tag_object)):
        die("local annotated seal tag does not bind the frozen seal")
    remote_raw, remote_rows = tag_remote_rows(
        repo, tag_name, environment, tool_records)
    expected_rows = {
        local_ref: tag_object,
        local_ref + "^{}": peeled,
    }
    if remote_rows and remote_rows != expected_rows:
        die("PERMANENT: GitHub seal tag name collides with another object")
    if not remote_rows:
        if now_utc_ms() >= seal_record["round_time_ms"]:
            die("holdout seal publication deadline passed before tag push")
        pushed = subprocess.run(
            (
                git, "-C", str(repo), "push", "origin",
                f"{local_ref}:{local_ref}",
            ),
            env=environment, text=True, capture_output=True,
            check=False, timeout=60,
        )
        if pushed.returncode:
            die(f"GitHub seal tag push failed: {pushed.stderr.strip()}")
        remote_raw, remote_rows = tag_remote_rows(
            repo, tag_name, environment, tool_records)
    confirmed_ms = now_utc_ms()
    if remote_rows != expected_rows:
        die("GitHub did not advertise the exact seal tag and peeled commit")
    return {
        "schema": "wirehair.wh2.holdout_seal_publication.v1",
        "status": (
            "CONFIRMED" if confirmed_ms < seal_record["round_time_ms"]
            else "ABANDONED_LATE"
        ),
        "phase": phase,
        "tag_name": tag_name,
        "tag_object": tag_object,
        "peeled_source_commit": peeled,
        "annotation_sha256": sha256_bytes(annotation.encode("utf-8")),
        "remote_rows_sha256": sha256_bytes(remote_raw),
        "remote_rows_ascii": remote_raw.decode("ascii"),
        "confirmed_ms": confirmed_ms,
    }


def r1_freeze_tag_annotation(
    prepare_record: dict[str, Any],
    prepare_sha256: str,
) -> str:
    return (
        "wirehair-wh2-r1-freeze-v1\n"
        f"prepare_sha256={prepare_sha256}\n"
        f"staged_seal_sha256={prepare_record['staged_seal_sha256']}\n"
        f"source_commit={prepare_record['source_commit']}\n"
        f"binary_sha256={prepare_record['binary_sha256']}\n"
        f"controller_sha256={prepare_record['script_sha256']}\n"
        f"contract_sha256={prepare_record['contract_sha256']}\n"
        f"protocol_sha256={prepare_record['protocol_sha256']}\n"
        f"route_context_sha256={prepare_record['route_context_sha256']}\n"
        f"development_seed_ledger_sha256="
        f"{prepare_record['development_seed_ledger_sha256']}\n"
        f"future_beacon_contract_sha256="
        f"{prepare_record['future_beacon_contract_sha256']}\n"
    )


def publish_r1_freeze_tag(
    repo: Path,
    prepare_record: dict[str, Any],
    prepare_sha256: str,
    tool_records: Any,
) -> dict[str, Any]:
    git = str(frozen_tool_path(tool_records, "git"))
    environment = newest_github_agent_environment(tool_records)
    tag_name = f"wh2-h12-kpreferred-r1-freeze-{prepare_sha256[:16]}"
    if not re.fullmatch(
            r"wh2-h12-kpreferred-r1-freeze-[0-9a-f]{16}", tag_name):
        die("internal R1 freeze tag name mismatch")
    annotation = r1_freeze_tag_annotation(prepare_record, prepare_sha256)
    local_ref = f"refs/tags/{tag_name}"
    local_exists = subprocess.run(
        (git, "-C", str(repo), "show-ref", "--verify", "--quiet", local_ref),
        check=False,
    ).returncode == 0
    if not local_exists:
        subprocess.run(
            (
                git, "-C", str(repo), "tag", "-a", tag_name,
                prepare_record["source_commit"], "-m", annotation,
            ),
            check=True, timeout=30,
        )
    tag_object = subprocess.check_output(
        (git, "-C", str(repo), "rev-parse", f"{tag_name}^{{tag}}"),
        text=True,
    ).strip()
    peeled = subprocess.check_output(
        (git, "-C", str(repo), "rev-parse", f"{tag_name}^{{}}"),
        text=True,
    ).strip()
    tag_bytes = subprocess.check_output(
        (git, "-C", str(repo), "cat-file", "tag", tag_object))
    try:
        tag_message = tag_bytes.split(b"\n\n", 1)[1].decode("utf-8")
    except (IndexError, UnicodeDecodeError):
        die("local annotated R1 freeze tag is malformed")
    if (peeled != prepare_record["source_commit"] or
            tag_message != annotation or
            not re.fullmatch(r"[0-9a-f]{40}", tag_object)):
        die("local R1 freeze tag does not bind the prepared campaign")
    remote_raw, remote_rows = tag_remote_rows(
        repo, tag_name, environment, tool_records)
    expected_rows = {local_ref: tag_object, local_ref + "^{}": peeled}
    if remote_rows and remote_rows != expected_rows:
        die("PERMANENT: R1 freeze tag name collides with another object")
    if not remote_rows:
        pushed = subprocess.run(
            (
                git, "-C", str(repo), "push", "origin",
                f"{local_ref}:{local_ref}",
            ),
            env=environment, text=True, capture_output=True,
            check=False, timeout=60,
        )
        if pushed.returncode:
            die(f"R1 freeze tag push failed: {pushed.stderr.strip()}")
        remote_raw, remote_rows = tag_remote_rows(
            repo, tag_name, environment, tool_records)
    if remote_rows != expected_rows:
        die("remote did not advertise the exact R1 freeze tag")
    return {
        "schema": "wirehair.wh2.r1_freeze_publication.v1",
        "status": "CONFIRMED",
        "prepare_sha256": prepare_sha256,
        "tag_name": tag_name,
        "tag_object": tag_object,
        "peeled_source_commit": peeled,
        "annotation_sha256": sha256_bytes(annotation.encode("utf-8")),
        "remote_rows_sha256": sha256_bytes(remote_raw),
        "remote_rows_ascii": remote_raw.decode("ascii"),
        "confirmed_ms": now_utc_ms(),
    }


def verify_r1_freeze_publication(
    result_dir: Path,
    prepare_record: dict[str, Any],
    contract: dict[str, Any],
) -> dict[str, Any]:
    prepare_path = result_dir / "prepare.json"
    prepare_bytes = common.stable_bytes(prepare_path)
    if prepare_bytes != common.json_bytes(prepare_record):
        die("PERMANENT: prepare record is not canonical JSON")
    prepare_sha256 = sha256_bytes(prepare_bytes)
    expected_tag = f"wh2-h12-kpreferred-r1-freeze-{prepare_sha256[:16]}"
    publication = load_fixed_json(
        result_dir / "frozen/r1_freeze_publication.json",
        "R1 freeze publication",
    )
    remote_ascii = publication.get("remote_rows_ascii")
    try:
        remote_sha256 = sha256_bytes(remote_ascii.encode("ascii")) \
            if isinstance(remote_ascii, str) else None
    except UnicodeEncodeError:
        remote_sha256 = None
    expected_keys = {
        "schema", "status", "prepare_sha256", "tag_name", "tag_object",
        "peeled_source_commit", "annotation_sha256", "remote_rows_sha256",
        "remote_rows_ascii", "confirmed_ms",
    }
    confirmed_ms = publication.get("confirmed_ms")
    prepared_ms = prepare_record.get("prepared_utc_ns", 0) // 1_000_000 \
        if isinstance(prepare_record.get("prepared_utc_ns"), int) else 0
    if (set(publication) != expected_keys or
            publication.get("schema") !=
                "wirehair.wh2.r1_freeze_publication.v1" or
            publication.get("status") != "CONFIRMED" or
            publication.get("prepare_sha256") != prepare_sha256 or
            publication.get("tag_name") != expected_tag or
            publication.get("peeled_source_commit") !=
                contract.get("source_commit") or
            not re.fullmatch(
                r"[0-9a-f]{40}", str(publication.get("tag_object", ""))) or
            not isinstance(confirmed_ms, int) or isinstance(confirmed_ms, bool) or
            confirmed_ms < prepared_ms or
            publication.get("remote_rows_sha256") != remote_sha256):
        die("PERMANENT: R1 freeze publication receipt mismatch")
    annotation = r1_freeze_tag_annotation(prepare_record, prepare_sha256)
    if publication.get("annotation_sha256") != \
            sha256_bytes(annotation.encode("utf-8")):
        die("PERMANENT: R1 freeze annotation binding changed")
    repo = Path(str(contract["source_repo"])).resolve(strict=True)
    tool_records = contract.get("system_tools")
    git = str(frozen_tool_path(tool_records, "git"))
    current_remote = subprocess.check_output(
        (git, "-C", str(repo), "remote", "get-url", "origin"), text=True,
    ).strip()
    if prepare_record.get("github_remote") != current_remote:
        die("PERMANENT: R1 freeze GitHub remote binding changed")
    tag_bytes = subprocess.check_output(
        (git, "-C", str(repo), "cat-file", "tag", publication["tag_object"]))
    try:
        tag_message = tag_bytes.split(b"\n\n", 1)[1].decode("utf-8")
    except (IndexError, UnicodeDecodeError):
        die("PERMANENT: local R1 freeze tag object is malformed")
    if tag_message != annotation:
        die("PERMANENT: local R1 freeze tag annotation changed")
    environment = newest_github_agent_environment(tool_records)
    raw, rows = tag_remote_rows(
        repo, expected_tag, environment, tool_records)
    expected_rows = {
        f"refs/tags/{expected_tag}": publication["tag_object"],
        f"refs/tags/{expected_tag}^{{}}": contract["source_commit"],
    }
    if rows != expected_rows:
        die("PERMANENT: confirmed remote R1 freeze tag changed or vanished")
    return {
        **publication,
        "reverified_remote_rows_sha256": sha256_bytes(raw),
        "reverified_ms": now_utc_ms(),
    }


def atomic_state_json(path: Path, value: dict[str, Any]) -> None:
    common.atomic_write(path, fixed_json_bytes(value))
    fsync_directory(path.parent)


def load_fixed_json(path: Path, context: str) -> dict[str, Any]:
    try:
        raw = common.stable_bytes(path)
        record = json.loads(raw)
    except (json.JSONDecodeError, UnicodeDecodeError, OSError):
        die(f"{context} is not valid JSON")
    if not isinstance(record, dict):
        die(f"{context} must be a JSON object")
    if fixed_json_bytes(record) != raw:
        die(f"{context} is not fixed-order canonical JSON")
    return record


def write_seal_attempt(
    result_dir: Path,
    spec: HoldoutPhaseSpec,
    manifest: bytes,
    seal_record: dict[str, Any],
    seal_bytes: bytes,
) -> tuple[Path, str]:
    seal_sha256 = sha256_bytes(seal_bytes)
    seal_parent = result_dir / "seals" / spec.seal_tree
    seal_parent.mkdir(parents=True, exist_ok=True)
    if seal_parent.is_symlink():
        die("seal attempt parent must not be a symlink")
    final = seal_parent / seal_sha256
    if final.exists():
        die("holdout seal digest already exists and may not be reused")
    staging = seal_parent / f".{seal_sha256}.partial-{os.getpid()}"
    staging.mkdir(mode=0o700, exist_ok=False)
    try:
        common.atomic_write(staging / "manifest.txt", manifest)
        common.atomic_write(staging / "seal.json", seal_bytes)
        common.atomic_write(
            staging / "seal.sha256", f"{seal_sha256}  seal.json\n".encode("ascii"))
        common.atomic_write(
            staging / "manifest.sha256",
            f"{sha256_bytes(manifest)}  manifest.txt\n".encode("ascii"),
        )
        if load_fixed_json(staging / "seal.json", "staged holdout seal") != seal_record:
            die("staged holdout seal changed")
        fsync_tree(staging)
        fsync_directory(seal_parent)
        os.rename(staging, final)
        fsync_directory(seal_parent)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return final, seal_sha256


def seal_state_path(result_dir: Path, phase: str) -> Path:
    state_dir = result_dir / "beacon" / phase
    state_dir.mkdir(mode=0o700, exist_ok=True)
    if state_dir.is_symlink():
        die("holdout state directory must not be a symlink")
    return state_dir / "state.json"


def load_seal_attempt(
    result_dir: Path,
    spec: HoldoutPhaseSpec,
    state: dict[str, Any],
) -> tuple[Path, bytes, dict[str, Any], str]:
    relative = state.get("seal_dir")
    seal_sha256 = state.get("seal_record_sha256")
    if (not isinstance(relative, str) or
            not re.fullmatch(
                rf"seals/{re.escape(spec.seal_tree)}/[0-9a-f]{{64}}", relative) or
            not isinstance(seal_sha256, str) or
            not re.fullmatch(r"[0-9a-f]{64}", seal_sha256) or
            not relative.endswith("/" + seal_sha256)):
        die("holdout controller state has an invalid seal path")
    candidate = result_dir / relative
    if candidate.is_symlink():
        die("holdout seal attempt must not be a symlink")
    attempt = candidate.resolve(strict=True)
    expected_parent = (result_dir / "seals" / spec.seal_tree).resolve(strict=True)
    if attempt.parent != expected_parent:
        die("holdout seal attempt escaped its immutable directory")
    entries = {path.name: path for path in attempt.iterdir()}
    required_names = {
        "manifest.txt", "seal.json", "seal.sha256", "manifest.sha256",
    }
    if (not required_names.issubset(entries) or
            not set(entries).issubset(required_names | {"publication.json"}) or
            any(path.is_symlink() or not path.is_file()
                for path in entries.values())):
        die("PERMANENT: holdout seal directory schema changed")
    seal_bytes = common.stable_bytes(attempt / "seal.json")
    if sha256_bytes(seal_bytes) != seal_sha256:
        die("PERMANENT: holdout seal record hash changed")
    seal_record = load_fixed_json(attempt / "seal.json", "holdout seal record")
    manifest = common.stable_bytes(attempt / "manifest.txt")
    if sha256_bytes(manifest) != seal_record.get("manifest_sha256"):
        die("PERMANENT: holdout seal manifest hash changed")
    if (common.stable_bytes(attempt / "seal.sha256") !=
            f"{seal_sha256}  seal.json\n".encode("ascii") or
            common.stable_bytes(attempt / "manifest.sha256") !=
            f"{sha256_bytes(manifest)}  manifest.txt\n".encode("ascii")):
        die("PERMANENT: holdout seal sidecar hash changed")
    return attempt, manifest, seal_record, seal_sha256


def validate_loaded_seal(
    spec: HoldoutPhaseSpec,
    manifest: bytes,
    seal_record: dict[str, Any],
    contract: dict[str, Any],
) -> None:
    seal_ms = seal_record.get("seal_ms")
    if (not isinstance(seal_ms, int) or isinstance(seal_ms, bool)):
        die("holdout seal timestamp is invalid")
    latest_bound = seal_record.get("verified_latest")
    latest_keys = {
        "schema", "latest_round", "canonical_beacon_sha256",
        "consensus_origins", "raw_response_sha256", "wrapper_result_sha256",
        "verified_ms",
    }
    all_origins = [f"https://{host}" for host in DRAND_RELAYS]
    if not isinstance(latest_bound, dict):
        die("holdout seal verified-latest bound is missing")
    latest_round = latest_bound.get("latest_round")
    consensus = latest_bound.get("consensus_origins")
    raw_hashes = latest_bound.get("raw_response_sha256")
    if (set(latest_bound) != latest_keys or
            latest_bound.get("schema") !=
                "wirehair.wh2.verified_latest_bound.v1" or
            not isinstance(latest_round, int) or isinstance(latest_round, bool) or
            not 1 <= latest_round <= (1 << 53) - 41 or
            not isinstance(consensus, list) or len(consensus) < 3 or
            len(consensus) != len(set(consensus)) or
            consensus != [origin for origin in all_origins if origin in consensus] or
            not isinstance(raw_hashes, dict) or
            not set(consensus).issubset(raw_hashes) or
            not set(raw_hashes).issubset(all_origins) or
            any(not isinstance(value, str) or
                not re.fullmatch(r"[0-9a-f]{64}", value)
                for value in raw_hashes.values()) or
            any(not isinstance(latest_bound.get(key), str) or
                not re.fullmatch(r"[0-9a-f]{64}", latest_bound[key])
                for key in ("canonical_beacon_sha256", "wrapper_result_sha256")) or
            not isinstance(latest_bound.get("verified_ms"), int) or
            isinstance(latest_bound.get("verified_ms"), bool)):
        die("holdout seal verified-latest bound is invalid")
    abandoned_round = seal_record.get("abandoned_round_exclusive")
    if (not isinstance(abandoned_round, int) or
            isinstance(abandoned_round, bool) or abandoned_round < 0):
        die("holdout seal abandoned-round bound is invalid")
    local_formula_round = target_beacon_round(seal_ms)
    expected_round = max(
        local_formula_round, latest_round + 40, abandoned_round + 1)
    expected_time = (
        DRAND_GENESIS_TIME * 1000 +
        (expected_round - 1) * DRAND_PERIOD_SECONDS * 1000
    )
    expected_chain = {
        "public_key": DRAND_PUBLIC_KEY,
        "period": DRAND_PERIOD_SECONDS,
        "genesis_time": DRAND_GENESIS_TIME,
        "hash": DRAND_CHAIN_HASH,
        "groupHash": DRAND_GROUP_HASH,
        "schemeID": "bls-unchained-g1-rfc9380",
        "metadata": {"beaconID": "quicknet"},
    }
    expected_bindings = {
        "phase": spec.phase,
        "manifest_sha256": sha256_bytes(manifest),
        "source_commit": contract["source_commit"],
        "dirty_diff_sha256": EMPTY_SHA256,
        "protocol_sha256": contract["protocol_sha256"],
        "binary_sha256": contract["binary_sha256"],
        "drand_wrapper_sha256": contract["drand_wrapper_sha256"],
        "drand_client_bundle_sha256": contract["drand_client_bundle_sha256"],
        "drand_package_sha256": contract["drand_package_sha256"],
        "drand_package_lock_sha256": contract["drand_package_lock_sha256"],
        "node_version": DRAND_NODE_VERSION,
        "npm_version": DRAND_NPM_VERSION,
        "drand_client_sri": DRAND_CLIENT_INTEGRITY,
        "chain_info": expected_chain,
        "origins": [f"https://{host}" for host in DRAND_RELAYS],
        "not_before_ms": seal_ms + 120_000,
        "local_formula_round": local_formula_round,
        "verified_latest": latest_bound,
        "verified_latest_minimum_round": latest_round + 40,
        "abandoned_round_exclusive": abandoned_round,
        "round": expected_round,
        "round_time_ms": expected_time,
    }
    expected_keys = {
        "schema", "github_remote", "seal_ms", *expected_bindings.keys(),
    }
    if (set(seal_record) != expected_keys or
            seal_record.get("schema") != "wirehair.wh2.holdout_seal.v1" or
            any(seal_record.get(key) != value
                for key, value in expected_bindings.items()) or
            seal_record.get("github_remote") != github_remote_url(
                Path(str(contract["source_repo"])).resolve(strict=True),
                frozen_tool_path(contract.get("system_tools"), "git"))):
        die("PERMANENT: holdout seal contract binding changed")


def validate_state_binding(
    state: dict[str, Any],
    spec: HoldoutPhaseSpec,
    manifest: bytes,
    seal_record: dict[str, Any],
    seal_sha256: str,
) -> None:
    if (state.get("schema") != "wirehair.wh2.holdout_controller_state.v1" or
            state.get("phase") != spec.phase or
            state.get("seal_record_sha256") != seal_sha256 or
            state.get("manifest_sha256") != sha256_bytes(manifest) or
            state.get("round") != seal_record.get("round") or
            state.get("round_time_ms") != seal_record.get("round_time_ms")):
        die("PERMANENT: holdout controller state binding changed")


def validate_rooted_state_binding(
    result_dir: Path,
    spec: HoldoutPhaseSpec,
    state: dict[str, Any],
) -> None:
    root_path = result_dir / spec.root_file
    if (state.get("status") != "ROOTED" or
            state.get("root_file") != spec.root_file or
            root_path.is_symlink() or not root_path.is_file() or
            state.get("root_file_sha256") != common.sha256_file(root_path)):
        die("PERMANENT: ROOTED state binding changed")


def seal_holdout(args: argparse.Namespace) -> int:
    result_dir = args.result_dir.resolve(strict=True)
    spec = HOLDOUT_PHASES[args.phase]
    with holdout_controller_lock(result_dir, args.phase):
        prepare_record, contract = verify_frozen_controller_runtime(result_dir)
        if args.phase == "h2":
            h1_state_path = result_dir / "beacon/h1/state.json"
            h1_state = load_fixed_json(h1_state_path, "H1 controller state")
            if h1_state.get("status") != "ROOTED":
                die("H2 sealing requires a fully rooted H1 phase")
            validate_rooted_state_binding(
                result_dir, HOLDOUT_PHASES["h1"], h1_state)
            h1_attempt, h1_manifest, h1_seal, h1_seal_sha = load_seal_attempt(
                result_dir, HOLDOUT_PHASES["h1"], h1_state)
            validate_loaded_seal(
                HOLDOUT_PHASES["h1"], h1_manifest, h1_seal, contract)
            validate_state_binding(
                h1_state, HOLDOUT_PHASES["h1"], h1_manifest,
                h1_seal, h1_seal_sha)
            verify_manifest_bytes(
                result_dir, HOLDOUT_PHASES["h1"], h1_manifest)
            verify_seal_publication(
                result_dir, h1_state, h1_attempt, h1_seal, h1_seal_sha,
                contract)
            verify_rooted_record(
                result_dir, HOLDOUT_PHASES["h1"], contract,
                h1_seal_sha, h1_seal,
            )
        state_path = seal_state_path(result_dir, args.phase)
        create_new = not state_path.exists()
        abandoned_round = 0
        if state_path.exists():
            state = load_fixed_json(state_path, f"{args.phase} controller state")
            if state.get("phase") != args.phase:
                die("holdout controller state phase mismatch")
            if state.get("status") in ("WAITING_UNTIL_T", "WAITING_FOR_QUORUM", "ROOTED"):
                attempt, manifest, seal_record, seal_sha256 = load_seal_attempt(
                    result_dir, spec, state)
                validate_loaded_seal(spec, manifest, seal_record, contract)
                validate_state_binding(
                    state, spec, manifest, seal_record, seal_sha256)
                verify_manifest_bytes(result_dir, spec, manifest)
                repo, commit = clean_source_commit(contract)
                if commit != seal_record.get("source_commit"):
                    die("holdout seal source commit mismatch")
                verify_seal_publication(
                    result_dir, state, attempt, seal_record, seal_sha256,
                    contract)
                if state.get("status") == "ROOTED":
                    verify_rooted_record(
                        result_dir, spec, contract, seal_sha256, seal_record)
                    validate_rooted_state_binding(result_dir, spec, state)
                print(canonical_json(state))
                return 0
            if str(state.get("status", "")).startswith("ABANDONED_"):
                abandoned_sha = state.get("seal_record_sha256")
                if (not isinstance(abandoned_sha, str) or
                        not re.fullmatch(r"[0-9a-f]{64}", abandoned_sha)):
                    die("abandoned holdout state has an invalid seal hash")
                (_abandoned_attempt, abandoned_manifest,
                 abandoned_seal, abandoned_loaded_sha) = load_seal_attempt(
                    result_dir, spec, state)
                validate_loaded_seal(
                    spec, abandoned_manifest, abandoned_seal, contract)
                validate_state_binding(
                    state, spec, abandoned_manifest,
                    abandoned_seal, abandoned_loaded_sha)
                abandoned_round = state.get("round", 0)
                if (not isinstance(abandoned_round, int) or
                        isinstance(abandoned_round, bool) or abandoned_round < 1):
                    die("abandoned holdout state has an invalid target round")
                abandoned = state_path.parent / (
                    f"abandoned-{abandoned_sha}.json")
                if abandoned.exists():
                    if common.stable_bytes(abandoned) != common.stable_bytes(state_path):
                        die("abandoned holdout state archive collision")
                else:
                    atomic_state_json(abandoned, state)
                create_new = True
            elif state.get("status") != "SEALED":
                die(f"holdout seal is not resumable from {state.get('status')!r}")
            else:
                attempt, manifest, seal_record, seal_sha256 = load_seal_attempt(
                    result_dir, spec, state)
                validate_loaded_seal(spec, manifest, seal_record, contract)
                validate_state_binding(
                    state, spec, manifest, seal_record, seal_sha256)
                verify_manifest_bytes(result_dir, spec, manifest)
                repo, commit = clean_source_commit(contract)
                if commit != seal_record.get("source_commit"):
                    die("holdout seal source commit mismatch")
        if create_new:
            for relative in spec.forbidden_before_seal:
                if (result_dir / relative).exists():
                    die(f"future phase input exists before {args.phase} seal: {relative}")
            verify_phase_ready_for_seal(
                result_dir, spec, prepare_record, contract)
            repo, commit = clean_source_commit(contract)
            remote_url = github_remote_url(
                repo, frozen_tool_path(contract.get("system_tools"), "git"))
            manifest = phase_manifest_bytes(result_dir, spec)
            latest_bound = obtain_verified_latest_bound(contract, result_dir)
            seal_ms = now_utc_ms()
            prepared_ms = int(prepare_record.get("prepared_utc_ns", 0)) // 1_000_000
            if seal_ms < prepared_ms:
                die("PERMANENT: UTC clock rolled back before holdout seal")
            seal_record, seal_bytes = build_seal_record(
                args.phase, sha256_bytes(manifest), contract, remote_url,
                seal_ms, latest_bound, abandoned_round)
            attempt, seal_sha256 = write_seal_attempt(
                result_dir, spec, manifest, seal_record, seal_bytes)
            verify_manifest_bytes(result_dir, spec, manifest)
            state = {
                "schema": "wirehair.wh2.holdout_controller_state.v1",
                "phase": args.phase,
                "status": "SEALED",
                "seal_dir": attempt.relative_to(result_dir).as_posix(),
                "seal_record_sha256": seal_sha256,
                "manifest_sha256": sha256_bytes(manifest),
                "round": seal_record["round"],
                "round_time_ms": seal_record["round_time_ms"],
            }
            atomic_state_json(state_path, state)
        if now_utc_ms() >= seal_record["round_time_ms"]:
            state["status"] = "ABANDONED_UNPUBLISHED_LATE"
            atomic_state_json(state_path, state)
            die("holdout seal reached T without confirmed tag publication")
        publication = publish_seal_tag(
            repo, args.phase, seal_sha256, sha256_bytes(manifest), seal_record,
            contract.get("system_tools"))
        publication["receipt_write_started_ms"] = now_utc_ms()
        atomic_state_json(attempt / "publication.json", publication)
        receipt_persisted_ms = now_utc_ms()
        if (publication["status"] != "CONFIRMED" or
                receipt_persisted_ms >= seal_record["round_time_ms"]):
            publication["status"] = "ABANDONED_LATE"
            publication["receipt_persisted_ms"] = receipt_persisted_ms
            atomic_state_json(attempt / "publication.json", publication)
            state["status"] = "ABANDONED_PUBLICATION_LATE"
            atomic_state_json(state_path, state)
            die("holdout seal tag confirmation completed at or after T")
        state.update({
            "status": "WAITING_UNTIL_T",
            "tag_name": publication["tag_name"],
            "tag_object": publication["tag_object"],
            "publication_confirmed_ms": publication["confirmed_ms"],
        })
        atomic_state_json(state_path, state)
        print(canonical_json(state))
        return 0


def verify_seal_publication(
    result_dir: Path,
    state: dict[str, Any],
    attempt: Path,
    seal_record: dict[str, Any],
    seal_sha256: str,
    contract: dict[str, Any],
) -> dict[str, Any]:
    publication = load_fixed_json(
        attempt / "publication.json", "holdout seal publication")
    phase = state["phase"]
    expected_tag = f"wh2-h12-kpreferred-{phase}-seal-{seal_sha256[:16]}"
    recorded_remote = publication.get("remote_rows_ascii")
    try:
        recorded_remote_sha256 = sha256_bytes(recorded_remote.encode("ascii")) \
            if isinstance(recorded_remote, str) else None
    except UnicodeEncodeError:
        recorded_remote_sha256 = None
    expected_keys = {
        "schema", "status", "phase", "tag_name", "tag_object",
        "peeled_source_commit", "annotation_sha256", "remote_rows_sha256",
        "remote_rows_ascii", "confirmed_ms", "receipt_write_started_ms",
    }
    confirmed_ms = publication.get("confirmed_ms")
    receipt_started_ms = publication.get("receipt_write_started_ms")
    if (set(publication) != expected_keys or
            publication.get("schema") !=
                "wirehair.wh2.holdout_seal_publication.v1" or
            publication.get("status") != "CONFIRMED" or
            publication.get("phase") != phase or
            publication.get("tag_name") != expected_tag or
            state.get("tag_name") != expected_tag or
            publication.get("tag_object") != state.get("tag_object") or
            state.get("publication_confirmed_ms") != confirmed_ms or
            publication.get("peeled_source_commit") !=
                seal_record.get("source_commit") or
            not isinstance(confirmed_ms, int) or isinstance(confirmed_ms, bool) or
            not seal_record["seal_ms"] <= confirmed_ms <
                seal_record["round_time_ms"] or
            not isinstance(receipt_started_ms, int) or
            isinstance(receipt_started_ms, bool) or
            not confirmed_ms <= receipt_started_ms <
                seal_record["round_time_ms"] or
            publication.get("remote_rows_sha256") != recorded_remote_sha256):
        die("PERMANENT: holdout seal publication receipt mismatch")
    repo = Path(str(contract["source_repo"]))
    repo = repo.resolve(strict=True)
    annotation = seal_tag_annotation(
        phase, seal_sha256, seal_record["manifest_sha256"], seal_record)
    if publication.get("annotation_sha256") != sha256_bytes(annotation.encode("utf-8")):
        die("PERMANENT: holdout seal tag annotation binding changed")
    tool_records = contract.get("system_tools")
    git = str(frozen_tool_path(tool_records, "git"))
    tag_bytes = subprocess.check_output(
        (git, "-C", str(repo), "cat-file", "tag", publication["tag_object"]))
    try:
        tag_message = tag_bytes.split(b"\n\n", 1)[1].decode("utf-8")
    except (IndexError, UnicodeDecodeError):
        die("PERMANENT: local holdout seal tag object is malformed")
    if tag_message != annotation:
        die("PERMANENT: local holdout seal tag annotation changed")
    environment = newest_github_agent_environment(tool_records)
    raw, rows = tag_remote_rows(
        repo, expected_tag, environment, tool_records)
    expected_rows = {
        f"refs/tags/{expected_tag}": publication["tag_object"],
        f"refs/tags/{expected_tag}^{{}}": seal_record["source_commit"],
    }
    if rows != expected_rows:
        die("PERMANENT: confirmed GitHub holdout seal tag changed or vanished")
    return {
        **publication,
        "reverified_remote_rows_sha256": sha256_bytes(raw),
        "reverified_ms": now_utc_ms(),
    }


def offline_verify_beacon(
    contract: dict[str, Any],
    result_dir: Path,
    beacon: dict[str, Any],
) -> dict[str, Any]:
    canonical = canonical_beacon_bytes(beacon)
    result = subprocess.run(
        (
            str(contract["node"]["path"]),
            str(result_dir / "frozen/wh2_drand_verify.cjs"),
            str(result_dir / "frozen/drand-client.cjs"),
            "offline-verify",
        ),
        input=canonical, capture_output=True, check=False, timeout=30,
    )
    if result.returncode or result.stderr:
        die("PERMANENT: stored Quicknet beacon failed offline BLS verification")
    try:
        record = json.loads(result.stdout)
    except (json.JSONDecodeError, UnicodeDecodeError):
        die("offline Quicknet verifier returned malformed JSON")
    if (record.get("schema") !=
            "wirehair.wh2.drand_quicknet_offline_verified.v1" or
            record.get("beacon") != beacon or
            record.get("canonical_sha256") != sha256_bytes(canonical)):
        die("offline Quicknet verifier record mismatch")
    return record


def earlier_h1_roots(result_dir: Path) -> tuple[str, ...]:
    record = load_fixed_json(
        result_dir / HOLDOUT_PHASES["h1"].root_file, "H1 holdout roots")
    roots = record.get("roots")
    if (not isinstance(roots, list) or len(roots) != H1_ROOT_COUNT or
            any(not isinstance(root, str) or
                not re.fullmatch(r"0x[0-9a-f]{16}", root) for root in roots)):
        die("H1 root file has invalid phase-root coverage")
    return tuple(roots)


def verify_rooted_record(
    result_dir: Path,
    spec: HoldoutPhaseSpec,
    contract: dict[str, Any],
    seal_sha256: str,
    seal_record: dict[str, Any],
) -> dict[str, Any]:
    record = load_fixed_json(result_dir / spec.root_file, f"{spec.phase} root file")
    beacon = record.get("beacon")
    expected_tag = f"wh2-h12-kpreferred-{spec.phase}-seal-{seal_sha256[:16]}"
    publication = load_fixed_json(
        result_dir / "seals" / spec.seal_tree / seal_sha256 / "publication.json",
        "rooted holdout publication",
    )
    wave_number = record.get("wave")
    wave_path = (result_dir / "beacon" / spec.phase / "waves" /
                 f"wave-{wave_number:06d}.json") if isinstance(wave_number, int) else None
    consensus = record.get("consensus_origins")
    if (record.get("schema") != "wirehair.wh2.holdout_roots.v1" or
            record.get("phase") != spec.phase or
            record.get("status") != "ROOTED" or
            record.get("seal_record_sha256") != seal_sha256 or
            record.get("manifest_sha256") !=
                seal_record.get("manifest_sha256") or
            record.get("round") != seal_record.get("round") or
            record.get("round_time_ms") != seal_record.get("round_time_ms") or
            not isinstance(beacon, dict) or
            beacon.get("round") != seal_record.get("round") or
            record.get("canonical_beacon_sha256") !=
                sha256_bytes(canonical_beacon_bytes(beacon)) or
            record.get("master_root_sha256") !=
                holdout_master_root(seal_sha256, beacon) or
            record.get("tag_name") != expected_tag or
            record.get("tag_object") != publication.get("tag_object") or
            not isinstance(consensus, list) or
            any(not isinstance(origin, str) for origin in consensus) or
            len(set(consensus)) < 3 or
            not set(consensus).issubset(
                {f"https://{host}" for host in DRAND_RELAYS}) or
            not isinstance(record.get("observations"), list) or
            not isinstance(wave_number, int) or wave_number < 1 or
            wave_path is None or not wave_path.is_file() or wave_path.is_symlink() or
            common.sha256_file(wave_path) != record.get("wave_record_sha256") or
            record.get("rooted_ms", 0) < seal_record.get("round_time_ms", 1)):
        die("PERMANENT: stored holdout root binding mismatch")
    persisted_wave = load_fixed_json(wave_path, "rooted holdout wave")
    wave_result = persisted_wave.get("result")
    if (not isinstance(wave_result, dict) or
            wave_result.get("status") != "QUORUM" or
            wave_result.get("beacon") != beacon or
            wave_result.get("consensus_origins") != consensus or
            wave_result.get("observations") != record.get("observations") or
            persisted_wave.get("wrapper_returncode") != 0):
        die("PERMANENT: stored root differs from its verified quorum wave")
    earlier = earlier_h1_roots(result_dir) if spec.phase == "h2" else ()
    expected_roots = derive_beacon_roots(
        spec.phase, seal_sha256, beacon, earlier_holdout_roots=earlier)
    if record.get("roots") != list(expected_roots):
        die("PERMANENT: stored holdout phase roots changed")
    offline_verify_beacon(contract, result_dir, beacon)
    return record


def existing_wave_count(
    waves_dir: Path,
    phase: str,
    seal_sha256: str,
) -> int:
    entries = sorted(waves_dir.iterdir())
    if any(not path.is_file() or path.is_symlink() or
           not re.fullmatch(r"wave-[0-9]{6}\.json", path.name)
           for path in entries):
        die("PERMANENT: unexpected entry in holdout wave ledger")
    paths = entries
    for expected, path in enumerate(paths, 1):
        if path.name != f"wave-{expected:06d}.json" or path.is_symlink():
            die("PERMANENT: holdout wave ledger has a gap or noncanonical path")
        record = load_fixed_json(path, "holdout wave record")
        if (record.get("schema") != "wirehair.wh2.holdout_wave_record.v1" or
                record.get("phase") != phase or
                record.get("seal_record_sha256") != seal_sha256 or
                record.get("wave") != expected):
            die("PERMANENT: holdout wave ledger binding changed")
    return len(paths)


def acquire_holdout(args: argparse.Namespace) -> int:
    result_dir = args.result_dir.resolve(strict=True)
    spec = HOLDOUT_PHASES[args.phase]
    with holdout_controller_lock(result_dir, args.phase):
        _prepare_record, contract = verify_frozen_controller_runtime(result_dir)
        clean_source_commit(contract)
        state_path = seal_state_path(result_dir, args.phase)
        state = load_fixed_json(state_path, f"{args.phase} controller state")
        if state.get("phase") != args.phase or state.get("status") not in (
                "WAITING_UNTIL_T", "WAITING_FOR_QUORUM", "ROOTED"):
            die("holdout phase is not ready for beacon acquisition")
        attempt, manifest, seal_record, seal_sha256 = load_seal_attempt(
            result_dir, spec, state)
        validate_loaded_seal(spec, manifest, seal_record, contract)
        validate_state_binding(state, spec, manifest, seal_record, seal_sha256)
        verify_manifest_bytes(result_dir, spec, manifest)
        publication = verify_seal_publication(
            result_dir, state, attempt, seal_record, seal_sha256, contract)
        root_path = result_dir / spec.root_file
        if root_path.exists():
            rooted = verify_rooted_record(
                result_dir, spec, contract, seal_sha256, seal_record)
            root_file_sha256 = common.sha256_file(root_path)
            if state["status"] == "ROOTED":
                validate_rooted_state_binding(result_dir, spec, state)
            if state["status"] != "ROOTED":
                state["status"] = "ROOTED"
                state["root_file"] = spec.root_file
                state["root_file_sha256"] = root_file_sha256
                atomic_state_json(state_path, state)
            print(canonical_json(rooted))
            return 0
        if state["status"] == "ROOTED":
            die("PERMANENT: ROOTED state lost its atomic root file")
        current_ms = now_utc_ms()
        if current_ms < seal_record["seal_ms"]:
            die("PERMANENT: UTC clock rolled back after holdout seal")
        if current_ms < seal_record["round_time_ms"]:
            state["status"] = "WAITING_UNTIL_T"
            atomic_state_json(state_path, state)
            if args.no_wait:
                print(canonical_json(state))
                return 75
            time.sleep((seal_record["round_time_ms"] - current_ms) / 1000)
        state["status"] = "WAITING_FOR_QUORUM"
        atomic_state_json(state_path, state)
        waves_dir = result_dir / "beacon" / args.phase / "waves"
        waves_dir.mkdir(parents=True, exist_ok=True)
        wave_number = existing_wave_count(
            waves_dir, args.phase, seal_sha256)
        started = time.monotonic()
        delays = (1, 2, 4, 8, 16, 30)
        while True:
            verify_manifest_bytes(result_dir, spec, manifest)
            wave_number += 1
            command = (
                str(contract["node"]["path"]),
                str(result_dir / "frozen/wh2_drand_verify.cjs"),
                str(result_dir / "frozen/drand-client.cjs"),
                str(seal_record["round"]),
            )
            try:
                wave = subprocess.run(
                    command, capture_output=True, check=False, timeout=30)
            except subprocess.TimeoutExpired as exc:
                wave = subprocess.CompletedProcess(
                    command, 124, exc.stdout or b"", exc.stderr or b"")
            try:
                wave_record = json.loads(wave.stdout)
            except (json.JSONDecodeError, UnicodeDecodeError):
                wave_record = {
                    "schema": "wirehair.wh2.drand_quicknet_wave.v1",
                    "status": "PERMANENT_WRAPPER_FAILURE",
                    "stdout_sha256": sha256_bytes(wave.stdout),
                    "stderr_sha256": sha256_bytes(wave.stderr),
                    "returncode": wave.returncode,
                }
            persisted_wave = {
                "schema": "wirehair.wh2.holdout_wave_record.v1",
                "phase": args.phase,
                "seal_record_sha256": seal_sha256,
                "wave": wave_number,
                "completed_ms": now_utc_ms(),
                "wrapper_returncode": wave.returncode,
                "wrapper_stderr_sha256": sha256_bytes(wave.stderr),
                "result": wave_record,
            }
            wave_path = waves_dir / f"wave-{wave_number:06d}.json"
            if wave_path.exists():
                die("holdout wave record would overwrite an earlier wave")
            atomic_state_json(wave_path, persisted_wave)
            status = wave_record.get("status")
            if status == "PERMANENT_VERIFIED_DISAGREEMENT" or (
                    status not in ("QUORUM", "TEMPORARY_NO_QUORUM")):
                state["status"] = "BLOCKED_PERMANENT"
                state["blocked_wave"] = wave_number
                atomic_state_json(state_path, state)
                die(f"holdout beacon controller permanently blocked at wave {wave_number}")
            if status == "QUORUM":
                if wave.returncode:
                    die("QUORUM wrapper result had a failing process status")
                beacon = wave_record.get("beacon")
                consensus_origins = wave_record.get("consensus_origins")
                observations = wave_record.get("observations")
                if (not isinstance(beacon, dict) or
                        beacon.get("round") != seal_record["round"] or
                        not isinstance(consensus_origins, list) or
                        any(not isinstance(origin, str)
                            for origin in consensus_origins) or
                        len(set(consensus_origins)) < 3 or
                        not set(consensus_origins).issubset(
                            {f"https://{host}" for host in DRAND_RELAYS}) or
                        not isinstance(observations, list) or
                        len(observations) != len(DRAND_RELAYS) or
                        any(not isinstance(item, dict) or
                            not isinstance(item.get("origin"), str)
                            for item in observations) or
                        {item.get("origin") for item in observations
                         } !=
                            {f"https://{host}" for host in DRAND_RELAYS}):
                    die("QUORUM wrapper result has invalid coverage")
                offline_verify = offline_verify_beacon(
                    contract, result_dir, beacon)
                earlier = earlier_h1_roots(result_dir) if args.phase == "h2" else ()
                roots = derive_beacon_roots(
                    args.phase, seal_sha256, beacon,
                    earlier_holdout_roots=earlier,
                )
                rooted = {
                    "schema": "wirehair.wh2.holdout_roots.v1",
                    "phase": args.phase,
                    "status": "ROOTED",
                    "seal_record_sha256": seal_sha256,
                    "manifest_sha256": seal_record["manifest_sha256"],
                    "tag_name": publication["tag_name"],
                    "tag_object": publication["tag_object"],
                    "round": seal_record["round"],
                    "round_time_ms": seal_record["round_time_ms"],
                    "canonical_beacon_sha256":
                        offline_verify["canonical_sha256"],
                    "beacon": beacon,
                    "master_root_sha256":
                        holdout_master_root(seal_sha256, beacon),
                    "roots": list(roots),
                    "consensus_origins": consensus_origins,
                    "observations": observations,
                    "wave": wave_number,
                    "wave_record_sha256": common.sha256_file(wave_path),
                    "rooted_ms": now_utc_ms(),
                }
                root_path.parent.mkdir(parents=True, exist_ok=True)
                if root_path.exists():
                    die("holdout root file appeared before atomic publication")
                atomic_state_json(root_path, rooted)
                verify_manifest_bytes(result_dir, spec, manifest)
                verify_rooted_record(
                    result_dir, spec, contract, seal_sha256, seal_record)
                state["status"] = "ROOTED"
                state["root_file"] = spec.root_file
                state["root_file_sha256"] = common.sha256_file(root_path)
                atomic_state_json(state_path, state)
                print(canonical_json(rooted))
                return 0
            if wave.returncode != 2:
                die("temporary no-quorum wrapper result had an unexpected status")
            state["last_temporary_wave"] = wave_number
            atomic_state_json(state_path, state)
            if time.monotonic() - started >= args.max_wait_seconds:
                print(canonical_json(state))
                return 75
            delay = delays[min(wave_number - 1, len(delays) - 1)]
            time.sleep(min(delay, max(0.0, args.max_wait_seconds -
                                      (time.monotonic() - started))))


@contextmanager
def atomic_result_directory(
    requested: Path,
) -> Iterator[tuple[Path, Path]]:
    """Build a campaign in a private sibling and publish it once complete."""
    absolute = requested.absolute()
    absolute.parent.mkdir(parents=True, exist_ok=True)
    parent = absolute.parent.resolve(strict=True)
    final = parent / absolute.name
    if not final.name or final.exists():
        die(f"result directory already exists: {final}")
    lock = parent / f".{final.name}.prepare.lock"
    try:
        lock_fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        die(f"prepare lock already exists: {lock}")
    staging = parent / f".{final.name}.prepare-{os.getpid()}"
    try:
        os.write(lock_fd, f"pid={os.getpid()}\n".encode("ascii"))
        os.fsync(lock_fd)
        fsync_directory(parent)
        staging.mkdir(mode=0o700, exist_ok=False)
        try:
            yield staging, final
            if final.exists():
                die(f"result directory appeared during prepare: {final}")
            fsync_tree(staging)
            fsync_directory(parent)
            os.rename(staging, final)
            fsync_directory(parent)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
    finally:
        os.close(lock_fd)
        try:
            lock.unlink()
            fsync_directory(parent)
        except FileNotFoundError:
            pass


def elf_build_id(binary: Path, readelf: Path) -> str:
    result = subprocess.run(
        (str(readelf), "--notes", str(binary)), text=True,
        capture_output=True, check=False, timeout=30,
    )
    matches = re.findall(r"^[ \t]*Build ID: ([0-9a-f]+)[ \t]*$",
                         result.stdout, flags=re.MULTILINE)
    if result.returncode or len(matches) != 1:
        die("frozen benchmark must expose one lowercase ELF build ID")
    return matches[0]


def prepare(args: argparse.Namespace) -> int:
    if args.acknowledge != "POST_SELECTION_H12_PREFERRED_ATTEMPT_V1":
        die("prepare requires the exact post-selection acknowledgement")
    script = Path(__file__).resolve(strict=True)
    allk = script.with_name("wh2_rank_floor_two_anchor_allk.py")
    helper = script.with_name("wh2_rank_floor_two_anchor_screen.py")
    campaign_module = script.with_name("wh2_preferred_attempt_campaign.py")
    holdout_module = script.with_name("wh2_preferred_attempt_holdout.py")
    timing_module = script.with_name("wh2_preferred_attempt_timing.py")
    drand_wrapper = script.with_name("wh2_drand_verify.cjs")
    drand_package = script.parent / "drand-verifier/package.json"
    drand_lock = script.parent / "drand-verifier/package-lock.json"
    tracked_inputs = (
        script, allk, helper, campaign_module, holdout_module, timing_module,
        drand_wrapper, drand_package, drand_lock,
    )
    system_tools = system_tool_identities()
    python_runtime = python_runtime_identity()
    repo, head, upstream = tracked_clean_sources(tracked_inputs, system_tools)
    remote_url = github_remote_url(
        repo, frozen_tool_path(system_tools, "git"))
    protocol, ledgers = derive_plan(args.source_result)
    binary_input = args.binary.absolute()
    if binary_input.is_symlink():
        die("--binary must not be a symlink")
    binary = binary_input.resolve(strict=True)
    build_dir = binary.parent.parent.resolve(strict=True)
    cache = build_dir / "CMakeCache.txt"
    if (binary != (build_dir / "codec/wirehair_v2_bench").resolve(strict=True) or
            common.cmake_source_directory(cache) != repo):
        die("--binary is not this source tree's wirehair_v2_bench")
    if args.workers != 128 or not 1 <= args.build_workers <= 128:
        die("campaign requires 128 workers and 1..128 build workers")
    build_command = [
        str(frozen_tool_path(system_tools, "taskset")), "-c", "0-127",
        str(frozen_tool_path(system_tools, "cmake")), "--build", str(build_dir),
        "--target", "wirehair_v2_bench", "--clean-first", "--parallel",
        str(args.build_workers),
    ]
    subprocess.run(build_command, check=True)
    if tracked_clean_sources(
            tracked_inputs, system_tools)[1:] != (head, upstream):
        die("repository changed during clean benchmark rebuild")
    if binary_input.is_symlink():
        die("clean rebuild replaced the binary with a symlink")
    binary = binary_input.resolve(strict=True)
    live_binary_sha256 = common.sha256_file(binary)
    readelf = frozen_tool_path(system_tools, "readelf")
    live_build_id = elf_build_id(binary, readelf)
    live_identity = q0_control_identity(binary)
    if (common.sha256_file(binary) != live_binary_sha256 or
            elf_build_id(binary, readelf) != live_build_id):
        die("benchmark binary changed during q0 identity validation")
    drand = obtain_drand_verifier(
        drand_wrapper, drand_package, drand_lock)
    if tracked_clean_sources(
            tracked_inputs, system_tools)[1:] != (head, upstream):
        die("repository changed during drand verifier freeze")
    host_runtime = host_runtime_identity()
    thermal = args.thermal.resolve(strict=True)
    thermal_mark = thermal_start(thermal, require_zero_edac=False)
    with atomic_result_directory(args.result_dir) as (result_dir, final_dir):
        frozen = result_dir / "frozen"
        frozen.mkdir()
        copied = {
            "script": frozen / script.name,
            "allk": frozen / allk.name,
            "helper": frozen / helper.name,
            "campaign_module": frozen / campaign_module.name,
            "holdout_module": frozen / holdout_module.name,
            "timing_module": frozen / timing_module.name,
            "drand_wrapper": frozen / drand_wrapper.name,
            "binary": frozen / "wirehair_v2_bench",
            "cmake_cache": frozen / "CMakeCache.txt",
            "drand_package": frozen / "drand-package.json",
            "drand_lock": frozen / "drand-package-lock.json",
            "drand_bundle": frozen / "drand-client.cjs",
        }
        for source, target in (
            (script, copied["script"]), (allk, copied["allk"]),
            (helper, copied["helper"]),
            (campaign_module, copied["campaign_module"]),
            (holdout_module, copied["holdout_module"]),
            (timing_module, copied["timing_module"]),
            (drand_wrapper, copied["drand_wrapper"]),
            (drand_package, copied["drand_package"]),
            (drand_lock, copied["drand_lock"]),
            (binary, copied["binary"]),
            (cache, copied["cmake_cache"]),
        ):
            shutil.copyfile(source, target)
        common.atomic_write(copied["drand_bundle"], drand["bundle"])
        copied["script"].chmod(0o755)
        copied["drand_wrapper"].chmod(0o755)
        copied["binary"].chmod(0o755)
        if (common.sha256_file(copied["binary"]) != live_binary_sha256 or
                elf_build_id(copied["binary"], readelf) != live_build_id):
            die("frozen benchmark differs from the q0-tested binary")
        identity = q0_control_identity(copied["binary"])
        if identity != live_identity:
            die("frozen benchmark q0 identity differs after copy")
        frozen_selftest = subprocess.run(
            (
                drand["node"], str(copied["drand_wrapper"]),
                str(copied["drand_bundle"]), "offline-selftest",
            ),
            text=True, capture_output=True, check=False, timeout=90,
        )
        if frozen_selftest.returncode or frozen_selftest.stderr:
            die(
                "frozen drand known-round self-test failed: "
                f"{frozen_selftest.stderr.strip()}"
            )
        try:
            frozen_selftest_record = json.loads(frozen_selftest.stdout)
        except json.JSONDecodeError:
            die("frozen drand known-round self-test returned malformed JSON")
        if (frozen_selftest_record.get("schema") !=
                "wirehair.wh2.drand_quicknet_offline_kat.v1" or
                frozen_selftest_record.get("beacon") != DRAND_KNOWN_ROUND or
                frozen_selftest_record.get("canonical_sha256") !=
                DRAND_CANONICAL_KNOWN_SHA256):
            die("frozen drand offline self-test record mismatch")
        for name, data in ledgers.items():
            common.atomic_write(frozen / name, data)
        common.atomic_json(frozen / "protocol.json", protocol)
        common.atomic_json(frozen / "q0_identity.json", identity)
        protocol_sha256 = common.sha256_file(frozen / "protocol.json")
        route_search_ledger = {
            "rounds": protocol["development"]["rounds"],
            "development_roots": protocol["development_roots"],
            "widths": list(WIDTHS), "schedules": list(SCHEDULES),
            "losses": list(LOSSES), "candidate_batch": 32,
        }
        route_context = {
            "schema": "wirehair.wh2.h12_preferred_attempt.route_context.v1",
            "source_commit": head,
            "binary_sha256": live_binary_sha256,
            "binary_elf_build_id": live_build_id,
            "protocol_sha256": protocol_sha256,
            "cli_schema": "preferredattempt.v2",
            "manifest_schema": "preferredattempt-route-manifest.v1",
            "q0_profile": {
                "precode_contract": 4, "packet_row_contract": 4,
                "completion": "mixed", "mix_count": 2,
                "packet_peel_seed_xor": "0x0",
                "two_anchor_cutoff": CUTOFF, "two_anchor_phase": 0,
                "seed_block_bytes": "actual",
            },
            "cohort_sha256": COHORT_SHA256,
            "bins_sha256": BINS_SHA256,
            "allk_bins_sha256": ALLK_BINS_SHA256,
            "development_seed_ledger_sha256": DEV_SEED_LEDGER_SHA256,
            "future_beacon_contract_sha256": sha256_bytes(
                common.json_bytes(protocol["future_holdout_beacon"])),
            "search_ledger_sha256": sha256_bytes(
                common.json_bytes(route_search_ledger)),
            "route_generator_command": [
                "frozen/wirehair_v2_bench", "preferredattempt",
                "--mode", "route", "--N", "<sealed-bin-K-list>",
                "--bb-list", "<round-width-list>", "--preferred-map",
                "<sealed-logical-route-map>", "--route-context-sha256",
                "<this-file-sha256>",
            ],
            "caller_supplied_context_forbidden": True,
        }
        common.atomic_json(frozen / "route_context.json", route_context)
        route_context_sha256 = common.sha256_file(
            frozen / "route_context.json")
        contract = {
            "schema": SCHEMA, "source_commit": head,
            "upstream_commit": upstream,
            "source_repo": str(repo), "build_command": build_command,
            "binary_sha256": live_binary_sha256,
            "binary_elf_build_id": live_build_id,
            "cmake_cache_sha256": common.sha256_file(copied["cmake_cache"]),
            "script_sha256": common.sha256_file(copied["script"]),
            "allk_sha256": common.sha256_file(copied["allk"]),
            "helper_sha256": common.sha256_file(copied["helper"]),
            "campaign_module_sha256": common.sha256_file(
                copied["campaign_module"]),
            "holdout_module_sha256": common.sha256_file(
                copied["holdout_module"]),
            "timing_module_sha256": common.sha256_file(
                copied["timing_module"]),
            "drand_wrapper_sha256": common.sha256_file(
                copied["drand_wrapper"]),
            "drand_package_sha256": common.sha256_file(
                copied["drand_package"]),
            "drand_package_lock_sha256": common.sha256_file(
                copied["drand_lock"]),
            "drand_client_bundle_sha256": common.sha256_file(
                copied["drand_bundle"]),
            "node": {
                "path": drand["node"], "version": drand["node_version"],
                "sha256": drand["node_sha256"],
            },
            "npm": {
                "path": drand["npm"], "version": drand["npm_version"],
                "sha256": drand["npm_sha256"],
            },
            "python": python_runtime,
            "system_tools": system_tools,
            "host_runtime": host_runtime,
            "drand_offline_selftest": frozen_selftest_record,
            "protocol_sha256": protocol_sha256,
            "q0_identity_sha256": common.sha256_file(
                frozen / "q0_identity.json"),
            "route_context_sha256": route_context_sha256,
            "thermal": str(thermal), "workers": args.workers,
            "cpu_set": list(range(128)), "timeout_seconds": args.timeout,
            "expected_recovery_jobs_max": 34_560,
            "expected_route_manifest_jobs_max": 1_664,
            "expected_all_process_jobs_max": 36_224,
            "expected_dev_logical_cells_max": 3_461_400,
            "expected_h1_logical_cells": 276_912,
            "expected_h2_logical_cells": 4_607_928,
            "saturated_timing_speed_claim_valid": False,
            "thermal_fields": list(THERMAL_FIELDS),
            "thermal_baseline": {
                "dev": thermal_mark["dev"], "ino": thermal_mark["ino"],
                "offset": thermal_mark["offset"],
                "edac_ce": thermal_mark["edac_ce"],
                "edac_ue": thermal_mark["edac_ue"],
                "monotonic_s": thermal_mark["monotonic_s"],
                "max_temperature_c": thermal_mark["max_temperature_c"],
                "row_sha256": sha256_bytes(thermal_mark["baseline_row"]),
            },
            "thermal_policy": {
                "cpu_limit_c": 90.0, "dimm_limit_c": 90.0,
                "timing_cpu_limit_c": 85.0,
                "timing_dimm_limit_c": 90.0,
                "consecutive_samples": 3, "stale_seconds": 5.0,
                "min_cpu_busy_pct": 95.0,
                "edac_ce_delta": 0, "edac_ue_delta": 0,
            },
            "production_route_table_present": False,
            "production_profile_id_present": False,
            "holdout_root_files_present": False,
        }
        common.atomic_json(frozen / "contract.json", contract)
        contract_sha256 = common.sha256_file(frozen / "contract.json")
        staged_files = [
            path for path in frozen.iterdir() if path.name != "staged.sha256"
        ]
        common.atomic_write(
            frozen / "staged.sha256",
            common.sha_manifest(result_dir, staged_files),
        )
        common.verify_sha_manifest(result_dir, frozen / "staged.sha256")
        if any((result_dir / relative).exists() for relative in (
                "holdout/h1_roots.json", "holdout/h2_roots.json",
                "seals/development_table_seal.json", "seals/h1_table_seal.json")):
            die("future holdout or preceding table-seal file exists during prepare")
        if tracked_clean_sources(
                tracked_inputs, system_tools)[1:] != (head, upstream):
            die("repository changed before atomic campaign publication")
        prepare_record = {
            "schema": "wirehair.wh2.h12_preferred_attempt.prepare.v1",
            "source_commit": head, "binary_sha256": contract["binary_sha256"],
            "script_sha256": contract["script_sha256"],
            "allk_sha256": contract["allk_sha256"],
            "helper_sha256": contract["helper_sha256"],
            "campaign_module_sha256": contract["campaign_module_sha256"],
            "holdout_module_sha256": contract["holdout_module_sha256"],
            "timing_module_sha256": contract["timing_module_sha256"],
            "binary_elf_build_id": live_build_id,
            "contract_sha256": contract_sha256,
            "protocol_sha256": contract["protocol_sha256"],
            "route_context_sha256": route_context_sha256,
            "cohort_sha256": COHORT_SHA256, "bins_sha256": BINS_SHA256,
            "allk_bins_sha256": ALLK_BINS_SHA256,
            "development_seed_ledger_sha256": DEV_SEED_LEDGER_SHA256,
            "future_beacon_contract_sha256":
                route_context["future_beacon_contract_sha256"],
            "drand_package_sha256": contract["drand_package_sha256"],
            "drand_package_lock_sha256":
                contract["drand_package_lock_sha256"],
            "drand_wrapper_sha256": contract["drand_wrapper_sha256"],
            "drand_client_bundle_sha256": DRAND_CLIENT_BUNDLE_SHA256,
            "node": contract["node"], "npm": contract["npm"],
            "python": contract["python"],
            "system_tools": contract["system_tools"],
            "host_runtime": contract["host_runtime"],
            "drand_offline_selftest": frozen_selftest_record,
            "q0_identity": identity,
            "q0_identity_sha256": contract["q0_identity_sha256"],
            "staged_seal_sha256": common.sha256_file(
                frozen / "staged.sha256"),
            "github_remote": remote_url,
            "result_dir": str(final_dir), "launch_state": "CLEAN_UNLAUNCHED",
            "holdout_root_files_present": False,
            "prepared_utc_ns": time.time_ns(),
        }
        common.atomic_json(result_dir / "prepare.json", prepare_record)
        prepare_sha256 = common.sha256_file(result_dir / "prepare.json")
        publication = publish_r1_freeze_tag(
            repo, prepare_record, prepare_sha256, contract["system_tools"])
        atomic_state_json(
            frozen / "r1_freeze_publication.json", publication)
        verify_r1_freeze_publication(result_dir, prepare_record, contract)
    print(canonical_json(prepare_record))
    return 0


def frozen_development_rounds(campaign: Any) -> tuple[Any, ...]:
    rounds = tuple(campaign.DEFAULT_ROUNDS)
    observed = tuple(
        (
            value.name, value.input_max, value.retain_max,
            tuple(value.root_indexes), tuple(value.widths),
            value.require_control_improvement,
        )
        for value in rounds
    )
    expected = tuple(
        (
            value.name, value.input_count, value.retain_count,
            value.root_indexes, value.widths, value.name == "r4",
        )
        for value in ROUNDS
    )
    if observed != expected:
        die("PERMANENT: campaign development rounds differ from protocol")
    return rounds


def expected_development_execution_files(
    result_dir: Path,
    campaign: Any,
    ledgers: Sequence[Any],
) -> set[Path]:
    expected: set[Path] = set()
    for ledger in ledgers:
        ledger.record()
        ledger_base = result_dir / "ledgers" / (
            f"{ledger.phase}-{ledger.kind}.json")
        completion = result_dir / "completed" / (
            f"{ledger.phase}-{ledger.kind}.json")
        expected.update((
            ledger_base, ledger_base.with_suffix(".json.sha256"),
            completion, completion.with_suffix(".json.sha256"),
        ))
        if ledger.jobs:
            expected.add(
                result_dir / "thermal" /
                f"{ledger.phase}-{ledger.kind}.csv")
        for job in ledger.jobs:
            expected.update(campaign.job_paths(result_dir, job).values())
    return {path.resolve(strict=True) for path in expected}


def exact_tree_files(result_dir: Path, top_names: Sequence[str]) -> set[Path]:
    files: set[Path] = set()
    for name in top_names:
        top = result_dir / name
        if top.is_symlink() or not top.is_dir():
            die(f"completed campaign lacks regular {name} directory")
        for path in top.rglob("*"):
            if path.is_symlink():
                die(f"completed campaign contains a symlink: {path}")
            if path.is_file():
                files.add(path.resolve(strict=True))
            elif not path.is_dir():
                die(f"completed campaign contains a nonregular path: {path}")
    return files


def write_exact_manifest_once(
    result_dir: Path,
    target: Path,
    files: Sequence[Path] | set[Path],
    campaign: Any,
) -> str:
    expected = {path.resolve(strict=True) for path in files}
    if not expected:
        die(f"cannot seal an empty manifest: {target}")
    encoded = common.sha_manifest(result_dir, expected)
    campaign.durable_write_once_or_same(target, encoded)
    verified = common.verify_sha_manifest(result_dir, target)
    if set(verified) != expected:
        die(f"manifest does not cover its exact declared file set: {target}")
    return common.sha256_file(target)


def run_development(args: argparse.Namespace) -> int:
    """Resume the frozen control/R1..R4 campaign and seal its selected table."""
    result_dir = args.result_dir.resolve(strict=True)
    prepare_record, contract = verify_frozen_controller_runtime(result_dir)
    campaign = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_campaign",
        "campaign_module_sha256")
    cohort, bins = load_frozen_cohort_bins(result_dir, campaign)
    rounds = frozen_development_rounds(campaign)
    protocol = load_json_object(
        result_dir / "frozen/protocol.json", "frozen protocol")
    expected_root_ledger = development_seed_ledger()
    if (common.sha256_file(result_dir / "frozen/protocol.json") !=
            contract.get("protocol_sha256") or
            protocol.get("development_roots") != expected_root_ledger):
        die("PERMANENT: frozen development root ledger changed")
    roots = tuple(expected_root_ledger["dev"])
    root_parent = prepare_record.get("staged_seal_sha256")
    if not isinstance(root_parent, str) or not re.fullmatch(
            r"[0-9a-f]{64}", root_parent):
        die("frozen development parent binding is malformed")
    route_context_sha256 = contract.get("route_context_sha256")
    if not isinstance(route_context_sha256, str) or not re.fullmatch(
            r"[0-9a-f]{64}", route_context_sha256):
        die("frozen route-context binding is malformed")
    binary = (result_dir / "frozen/wirehair_v2_bench").resolve(strict=True)
    thermal = Path(str(contract.get("thermal", ""))).resolve(strict=True)
    cpus = tuple(contract.get("cpu_set", ()))
    workers = contract.get("workers")
    timeout = contract.get("timeout_seconds")
    runner = campaign.JobRunner(
        result_dir, binary, cpus, workers, timeout)
    evidence = campaign.DevelopmentEvidence(SCHEDULES, LOSSES)
    ledgers: list[Any] = []

    control = campaign.build_control_ledger(
        bins, cohort, roots, root_parent, SCHEDULES, LOSSES, WIDTHS)
    ledgers.append(control)
    evidence.add_completed(runner.run_ledger(control, thermal))

    development = result_dir / "development"
    development.mkdir(mode=0o700, exist_ok=True)
    rounds_dir = development / "rounds"
    rounds_dir.mkdir(mode=0o700, exist_ok=True)
    previous_record: dict[str, Any] | None = None
    previous_sha256 = root_parent
    round_rows: list[dict[str, Any]] = []
    for spec in rounds:
        derived = campaign.derive_round_input(
            rounds, spec, cohort, root_parent, previous_record)
        if derived.parent_sha256 != previous_sha256:
            die(f"{spec.name} no-backfill parent binding changed")
        route_ledger = campaign.build_route_ledger(
            spec, bins, cohort, derived.survivors,
            derived.parent_sha256, route_context_sha256)
        ledgers.append(route_ledger)
        route_completed = runner.run_ledger(route_ledger, thermal)
        evidence.add_completed(route_completed)
        route_artifacts = campaign.route_artifacts_from_completed(
            result_dir, route_completed)
        candidate_ledger = campaign.build_candidate_ledger(
            spec, route_ledger, route_artifacts, roots,
            route_ledger.sha256(), SCHEDULES, LOSSES)
        ledgers.append(candidate_ledger)
        candidate_completed = runner.run_ledger(candidate_ledger, thermal)
        evidence.add_completed(candidate_completed)
        survivor_record = campaign.analyze_round(
            result_dir, evidence, rounds, spec, cohort,
            derived.survivors, derived.parent_sha256)
        survivor_path = rounds_dir / f"{spec.name}.survivors.json"
        survivor_sha256 = campaign.write_survivor_ledger(
            result_dir, survivor_path, survivor_record)
        if survivor_sha256 != campaign.sha256_bytes(
                campaign.canonical_json_bytes(survivor_record)):
            die(f"{spec.name} survivor ledger hash changed while publishing")
        retained = campaign.survivors_from_record(
            survivor_record, spec.name, cohort, spec.retain_max)
        round_rows.append({
            "round": spec.name,
            "survivor_sha256": survivor_sha256,
            "input_candidates": survivor_record["input_candidates"],
            "retained_candidates": survivor_record["retained_candidates"],
            "routed_K": sum(bool(values) for values in retained.values()),
            "logical_cells": survivor_record["current_round_logical_cells"],
            "physical_candidate_solves":
                survivor_record["current_round_physical_candidate_solves"],
            "alias_cells": survivor_record["current_round_alias_cells"],
        })
        previous_record = survivor_record
        previous_sha256 = survivor_sha256
    if previous_record is None:
        die("development campaign produced no round record")
    final_survivors = campaign.survivors_from_record(
        previous_record, "r4", cohort, 1)
    global_objective = campaign.global_table_objective_report(
        result_dir, evidence, rounds, cohort, previous_record)

    expected_execution = expected_development_execution_files(
        result_dir, campaign, ledgers)
    actual_execution = exact_tree_files(
        result_dir, ("jobs", "ledgers", "completed", "thermal"))
    if actual_execution != expected_execution:
        missing = len(expected_execution - actual_execution)
        extra = len(actual_execution - expected_execution)
        die(
            "development execution artifact set mismatch "
            f"(missing={missing}, extra={extra})")
    phase_sha256 = write_exact_manifest_once(
        result_dir, development / "phase_complete.sha256",
        expected_execution, campaign)

    global_record = campaign.sealed_record(
        "wirehair.wh2.h12_preferred_attempt.global_table_objective.v1",
        {
            "phase_seal_sha256": phase_sha256,
            "final_survivor_sha256": previous_sha256,
            "objective": global_objective,
        },
    )
    global_path = development / "global_table_objective.json"
    global_sha256 = campaign.write_hashed_artifact(
        global_path, campaign.canonical_json_bytes(global_record))

    route_rows = [
        {
            "K": K,
            "route_status": "preferred" if final_survivors[K] else "control",
            "preferred_attempt": (
                final_survivors[K][0] if final_survivors[K] else None),
        }
        for K in cohort
    ]
    table = campaign.sealed_record(
        "wirehair.wh2.h12_preferred_attempt.development_table.v1",
        {
            "source_commit": prepare_record["source_commit"],
            "staged_seal_sha256": root_parent,
            "protocol_sha256": contract["protocol_sha256"],
            "route_context_sha256": route_context_sha256,
            "phase_seal_sha256": phase_sha256,
            "final_survivor_sha256": previous_sha256,
            "global_table_objective_sha256": global_sha256,
            "global_table_objective": global_objective,
            "K_count": len(cohort),
            "routed_K": sum(bool(values) for values in final_survivors.values()),
            "control_K": sum(not values for values in final_survivors.values()),
            "routes": route_rows,
        },
    )
    table_path = development / "frozen_table.json"
    table_sha256 = campaign.write_hashed_artifact(
        table_path, campaign.canonical_json_bytes(table))
    holdout = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_holdout",
        "holdout_module_sha256")
    selected_table = {
        K: (final_survivors[K][0] if final_survivors[K] else None)
        for K in cohort
    }
    route_cache_bindings: dict[int, Any] = {}
    route_cache_rows: list[dict[str, Any]] = []
    for bin_index in range(BIN_COUNT):
        Ks = tuple(bins[bin_index])
        source_records: list[Any] = []
        for K in Ks:
            p = selected_table[K]
            for width in WIDTHS:
                try:
                    a0 = evidence.control_canonical[(K, width)]
                except KeyError:
                    die("development control canonical map is incomplete")
                if p is None:
                    record = campaign.RouteRecord(
                        K, width, "control", -1, a0, a0,
                        1, 0, 1, 0, a0 + 1, 0)
                else:
                    try:
                        old = evidence.routes[(K, width, p)]
                    except KeyError:
                        die("selected development route evidence is incomplete")
                    record = campaign.RouteRecord(
                        old.K, old.width, old.route_status,
                        old.preferred_attempt, old.canonical_attempt,
                        old.actual_attempt, old.preferred_valid,
                        old.fallback, old.no_op, old.direct,
                        a0 + 1, old.preferred_probe_solves)
                campaign.validate_route_record_value(record)
                source_records.append(record)
        relative = f"development/route_cache/bin-{bin_index:03d}.csv"
        binding = holdout.publish_derived_route_cache(
            result_dir, relative, source_records, Ks, WIDTHS,
            {K: selected_table[K] for K in Ks}, route_context_sha256)
        route_cache_bindings[bin_index] = binding
        route_cache_rows.append({
            "bin": bin_index, "K_count": len(Ks),
            "path": binding.relative_path, "sha256": binding.sha256,
        })
    if set(route_cache_bindings) != set(range(BIN_COUNT)):
        die("development route-cache index is incomplete")
    route_cache_index = campaign.sealed_record(
        "wirehair.wh2.h12_preferred_attempt.development_route_cache_index.v1",
        {
            "frozen_table_sha256": table_sha256,
            "route_context_sha256": route_context_sha256,
            "cache_count": len(route_cache_rows),
            "caches": route_cache_rows,
        },
    )
    route_cache_index_path = development / "route_cache/index.json"
    route_cache_index_sha256 = campaign.write_hashed_artifact(
        route_cache_index_path,
        campaign.canonical_json_bytes(route_cache_index))
    summary = campaign.sealed_record(
        "wirehair.wh2.h12_preferred_attempt.development_summary.v1",
        {
            "phase_seal_sha256": phase_sha256,
            "frozen_table_sha256": table_sha256,
            "route_cache_index_sha256": route_cache_index_sha256,
            "global_table_objective_sha256": global_sha256,
            "global_table_objective": global_objective,
            "rounds": round_rows,
            "routed_K": table["routed_K"],
            "control_K": table["control_K"],
            "raw_architecture": "H12 q0 two-anchor selected before seed tuning",
            "seed_fixes_used_for_architecture_selection": False,
        },
    )
    summary_path = development / "summary.json"
    campaign.write_hashed_artifact(
        summary_path, campaign.canonical_json_bytes(summary))
    analysis_files = {
        path.resolve(strict=True)
        for path in development.rglob("*")
        if path.is_file() and path.name not in (
            "phase_complete.sha256", "analysis_complete.sha256")
    }
    analysis_sha256 = write_exact_manifest_once(
        result_dir, development / "analysis_complete.sha256",
        analysis_files, campaign)
    complete = {
        "schema": "wirehair.wh2.h12_preferred_attempt.development_complete.v1",
        "result_dir": str(result_dir),
        "phase_seal_sha256": phase_sha256,
        "analysis_seal_sha256": analysis_sha256,
        "frozen_table_sha256": table_sha256,
        "routed_K": table["routed_K"],
        "control_K": table["control_K"],
    }
    print(canonical_json(complete))
    return 0


def verify_named_hash_sidecar(path: Path) -> str:
    digest = common.sha256_file(path.resolve(strict=True))
    expected = f"{digest}  {path.name}\n".encode("ascii")
    sidecar = path.with_suffix(path.suffix + ".sha256")
    if common.stable_bytes(sidecar) != expected:
        die(f"hashed artifact sidecar changed: {path}")
    return digest


def load_development_table_and_caches(
    result_dir: Path,
    campaign: Any,
    holdout: Any,
    cohort: Sequence[int],
    bins: dict[int, tuple[int, ...]],
    route_context_sha256: str,
) -> tuple[dict[int, int | None], dict[int, Any], str]:
    table_path = result_dir / "development/frozen_table.json"
    table_record = campaign.load_canonical_object(
        table_path, "development frozen table")
    campaign.verify_sealed_record(
        table_record,
        "wirehair.wh2.h12_preferred_attempt.development_table.v1")
    exact_fields = {
        "schema", "source_commit", "staged_seal_sha256", "protocol_sha256",
        "route_context_sha256", "phase_seal_sha256",
        "final_survivor_sha256", "global_table_objective_sha256",
        "global_table_objective", "K_count", "routed_K", "control_K",
        "routes", "self_sha256_excluding_field",
    }
    rows = table_record.get("routes")
    if (set(table_record) != exact_fields or
            table_record.get("route_context_sha256") != route_context_sha256 or
            not isinstance(rows, list)):
        die("development frozen-table schema/binding changed")
    objective = campaign.validate_global_table_objective_report(
        table_record.get("global_table_objective"))
    global_path = result_dir / "development/global_table_objective.json"
    global_record = campaign.load_canonical_object(
        global_path, "development global table objective")
    campaign.verify_sealed_record(
        global_record,
        "wirehair.wh2.h12_preferred_attempt.global_table_objective.v1")
    global_sha256 = verify_named_hash_sidecar(global_path)
    if (set(global_record) != {
            "schema", "phase_seal_sha256", "final_survivor_sha256",
            "objective", "self_sha256_excluding_field"} or
            global_record.get("phase_seal_sha256") !=
                table_record.get("phase_seal_sha256") or
            global_record.get("final_survivor_sha256") !=
                table_record.get("final_survivor_sha256") or
            global_record.get("objective") != objective or
            table_record.get("global_table_objective_sha256") != global_sha256):
        die("development global-table objective binding changed")
    table: dict[int, int | None] = {}
    for row in rows:
        if not isinstance(row, dict) or set(row) != {
                "K", "route_status", "preferred_attempt"}:
            die("development frozen-table row schema changed")
        K = row.get("K")
        status = row.get("route_status")
        p = row.get("preferred_attempt")
        if (not isinstance(K, int) or isinstance(K, bool) or K in table or
                status not in ("preferred", "control") or
                ((status == "control") != (p is None)) or
                (p is not None and (
                    not isinstance(p, int) or isinstance(p, bool) or
                    not 0 <= p <= 255))):
            die("development frozen-table row is malformed")
        table[K] = p
    if (tuple(table) != tuple(cohort) or
            table_record.get("K_count") != len(cohort) or
            table_record.get("routed_K") !=
                sum(value is not None for value in table.values()) or
            table_record.get("control_K") !=
                sum(value is None for value in table.values())):
        die("development frozen-table coverage/arithmetic changed")
    table_sha256 = verify_named_hash_sidecar(table_path)
    index_path = result_dir / "development/route_cache/index.json"
    index = campaign.load_canonical_object(
        index_path, "development route-cache index")
    campaign.verify_sealed_record(
        index,
        "wirehair.wh2.h12_preferred_attempt.development_route_cache_index.v1")
    if (set(index) != {
            "schema", "frozen_table_sha256", "route_context_sha256",
            "cache_count", "caches", "self_sha256_excluding_field"} or
            index.get("frozen_table_sha256") != table_sha256 or
            index.get("route_context_sha256") != route_context_sha256 or
            index.get("cache_count") != BIN_COUNT or
            not isinstance(index.get("caches"), list)):
        die("development route-cache index binding changed")
    bindings: dict[int, Any] = {}
    for expected_bin, row in enumerate(index["caches"]):
        expected_path = f"development/route_cache/bin-{expected_bin:03d}.csv"
        if (not isinstance(row, dict) or set(row) != {
                "bin", "K_count", "path", "sha256"} or
                row.get("bin") != expected_bin or
                row.get("K_count") != len(bins[expected_bin]) or
                row.get("path") != expected_path):
            die("development route-cache index row changed")
        cache_path = result_dir / expected_path
        digest = verify_named_hash_sidecar(cache_path)
        if row.get("sha256") != digest:
            die("development route-cache digest changed")
        binding = holdout.RouteCacheBinding(expected_path, digest)
        binding.validate()
        holdout.parse_route_manifest(
            common.stable_bytes(cache_path), bins[expected_bin], WIDTHS,
            {K: table[K] for K in bins[expected_bin]},
            route_context_sha256, True)
        bindings[expected_bin] = binding
    verify_named_hash_sidecar(index_path)
    return table, bindings, table_sha256


def verify_stem_hash_sidecar(path: Path) -> str:
    """Verify the holdout module's ``stem.sha256`` sidecar convention."""
    digest = common.sha256_file(path.resolve(strict=True))
    expected = f"{digest}  {path.name}\n".encode("ascii")
    if common.stable_bytes(path.with_suffix(".sha256")) != expected:
        die(f"named record sidecar changed: {path}")
    return digest


def load_h1_table_and_caches(
    result_dir: Path,
    campaign: Any,
    holdout: Any,
    cohort: Sequence[int],
    bins: dict[int, tuple[int, ...]],
    development_table: dict[int, int | None],
    development_table_sha256: str,
    route_context_sha256: str,
) -> tuple[dict[int, int | None], dict[int, Any], str, str]:
    """Load the H1-pruned table and its exact four-width cache inventory."""
    h1_table_path = result_dir / "h1/frozen_table.json"
    h1_record = campaign.load_canonical_object(
        h1_table_path, "H1 frozen table")
    h1_table = holdout.parse_h1_table(h1_record, strict_domain=True)
    if tuple(h1_table) != tuple(cohort):
        die("H1 frozen table differs from the exact active cohort")
    for K in cohort:
        if (h1_table[K] is not None and
                h1_table[K] != development_table[K]):
            die("H1 table substituted a route instead of pruning to control")
    h1_table_sha256 = common.sha256_file(h1_table_path.resolve(strict=True))

    phase_path = result_dir / "h1/phase_complete.json"
    phase_sha256 = verify_stem_hash_sidecar(phase_path)
    analysis_path = result_dir / "h1/analysis_complete.json"
    analysis = campaign.load_canonical_object(
        analysis_path, "H1 analysis completion")
    campaign.verify_sealed_record(
        analysis,
        "wirehair.wh2.h12_preferred_attempt.holdout.h1_analysis.v1")
    analysis_sha256 = verify_stem_hash_sidecar(analysis_path)
    evidence = analysis.get("evidence_binding")
    expected_evidence = {
        "phase_completion_sha256": phase_sha256,
        "root_file_sha256": common.sha256_file(
            result_dir / HOLDOUT_PHASES["h1"].root_file),
        "parent_table_sha256": development_table_sha256,
        "route_context_sha256": route_context_sha256,
        "input_table_sha256": holdout.table_semantic_sha256(
            cohort, development_table),
    }
    if (analysis.get("accepted") is not True or
            analysis.get("gates") != {"errors_zero": True} or
            analysis.get("frozen_table") != h1_record or
            analysis.get("frozen_table_file_sha256") != h1_table_sha256 or
            analysis.get("input_table_sha256") !=
                expected_evidence["input_table_sha256"] or
            analysis.get("output_table_sha256") !=
                holdout.table_semantic_sha256(cohort, h1_table) or
            evidence != expected_evidence):
        die("H1 analysis/table evidence binding changed")

    index_path = result_dir / "h1/route_cache/index.json"
    index = campaign.load_canonical_object(
        index_path, "H1 route-cache index")
    campaign.verify_sealed_record(
        index,
        "wirehair.wh2.h12_preferred_attempt.h1_route_cache_index.v1")
    if (set(index) != {
            "schema", "frozen_table_sha256", "route_context_sha256",
            "cache_count", "caches", "self_sha256_excluding_field"} or
            index.get("frozen_table_sha256") != h1_table_sha256 or
            index.get("route_context_sha256") != route_context_sha256 or
            index.get("cache_count") != BIN_COUNT or
            not isinstance(index.get("caches"), list)):
        die("H1 route-cache index binding changed")
    bindings: dict[int, Any] = {}
    for expected_bin, row in enumerate(index["caches"]):
        expected_path = f"h1/route_cache/bin-{expected_bin:03d}.csv"
        if (not isinstance(row, dict) or set(row) != {
                "bin", "K_count", "path", "sha256"} or
                row.get("bin") != expected_bin or
                row.get("K_count") != len(bins[expected_bin]) or
                row.get("path") != expected_path):
            die("H1 route-cache index row changed")
        cache_path = result_dir / expected_path
        digest = verify_named_hash_sidecar(cache_path)
        if row.get("sha256") != digest:
            die("H1 route-cache digest changed")
        binding = holdout.RouteCacheBinding(expected_path, digest)
        binding.validate()
        holdout.parse_route_manifest(
            common.stable_bytes(cache_path), bins[expected_bin], WIDTHS,
            {K: h1_table[K] for K in bins[expected_bin]},
            route_context_sha256, True)
        bindings[expected_bin] = binding
    verify_named_hash_sidecar(index_path)
    return h1_table, bindings, h1_table_sha256, analysis_sha256


def verify_exact_sha_manifest_coverage(
    result_dir: Path,
    manifest_path: Path,
    expected_files: set[Path],
    context: str,
) -> str:
    """Require one canonical SHA manifest to cover exactly a frozen inventory."""
    expected = {path.resolve(strict=True) for path in expected_files}
    encoded = common.sha_manifest(result_dir, expected)
    if common.stable_bytes(manifest_path) != encoded:
        die(f"{context} is not the canonical exact file inventory")
    verified = common.verify_sha_manifest(result_dir, manifest_path)
    if set(verified) != expected:
        die(f"{context} file coverage changed")
    return common.sha256_file(manifest_path)


def verify_development_seal_readiness(
    result_dir: Path,
    prepare_record: dict[str, Any],
    contract: dict[str, Any],
) -> None:
    """Validate the completed development table before revealing H1."""
    campaign = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_campaign",
        "campaign_module_sha256")
    holdout = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_holdout",
        "holdout_module_sha256")
    cohort, bins = load_frozen_cohort_bins(result_dir, campaign)
    route_context_sha256 = contract.get("route_context_sha256")
    if (not isinstance(route_context_sha256, str) or
            not re.fullmatch(r"[0-9a-f]{64}", route_context_sha256)):
        die("frozen route-context binding is malformed")
    _table, _caches, _table_sha256 = load_development_table_and_caches(
        result_dir, campaign, holdout, cohort, bins, route_context_sha256)

    phase_path = result_dir / "development/phase_complete.sha256"
    phase_inventory = exact_tree_files(
        result_dir, ("jobs", "ledgers", "completed", "thermal"))
    phase_sha256 = verify_exact_sha_manifest_coverage(
        result_dir, phase_path, phase_inventory,
        "development phase-completion manifest")

    analysis_path = result_dir / "development/analysis_complete.sha256"
    analysis_inventory = exact_tree_files(result_dir, ("development",))
    analysis_inventory.remove(phase_path.resolve(strict=True))
    analysis_inventory.remove(analysis_path.resolve(strict=True))
    verify_exact_sha_manifest_coverage(
        result_dir, analysis_path, analysis_inventory,
        "development analysis-completion manifest")

    table_record = campaign.load_canonical_object(
        result_dir / "development/frozen_table.json",
        "development frozen table")
    expected_bindings = {
        "source_commit": contract.get("source_commit"),
        "staged_seal_sha256": prepare_record.get("staged_seal_sha256"),
        "protocol_sha256": contract.get("protocol_sha256"),
        "route_context_sha256": route_context_sha256,
        "phase_seal_sha256": phase_sha256,
    }
    if any(table_record.get(key) != value
           for key, value in expected_bindings.items()):
        die("development table does not bind the exact frozen campaign")


def h1_controller_complete_record(
    campaign: Any,
    *,
    phase_completion_sha256: str,
    analysis_complete_sha256: str,
    root_file_sha256: str,
    parent_table_sha256: str,
    frozen_table_sha256: str,
    route_cache_index_sha256: str,
    route_context_sha256: str,
    input_routed_K: int,
    retained_routed_K: int,
) -> dict[str, Any]:
    return campaign.sealed_record(
        "wirehair.wh2.h12_preferred_attempt.h1_controller_complete.v1",
        {
            "phase_completion_sha256": phase_completion_sha256,
            "analysis_complete_sha256": analysis_complete_sha256,
            "root_file_sha256": root_file_sha256,
            "parent_table_sha256": parent_table_sha256,
            "frozen_table_sha256": frozen_table_sha256,
            "route_cache_index_sha256": route_cache_index_sha256,
            "route_context_sha256": route_context_sha256,
            "input_routed_K": input_routed_K,
            "retained_routed_K": retained_routed_K,
            "cache_count": BIN_COUNT,
        },
    )


def verify_h1_seal_readiness(
    result_dir: Path,
    contract: dict[str, Any],
) -> None:
    """Validate the terminal H1 table/cache/evidence state before H2."""
    campaign = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_campaign",
        "campaign_module_sha256")
    holdout = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_holdout",
        "holdout_module_sha256")
    cohort, bins = load_frozen_cohort_bins(result_dir, campaign)
    route_context_sha256 = contract.get("route_context_sha256")
    if (not isinstance(route_context_sha256, str) or
            not re.fullmatch(r"[0-9a-f]{64}", route_context_sha256)):
        die("frozen route-context binding is malformed")
    development_table, development_caches, development_table_sha256 = \
        load_development_table_and_caches(
            result_dir, campaign, holdout, cohort, bins,
            route_context_sha256)
    h1_table, h1_caches, h1_table_sha256, analysis_sha256 = \
        load_h1_table_and_caches(
            result_dir, campaign, holdout, cohort, bins,
            development_table, development_table_sha256,
            route_context_sha256)

    roots = earlier_h1_roots(result_dir)
    root_file_sha256 = common.sha256_file(
        result_dir / HOLDOUT_PHASES["h1"].root_file)
    ledger = holdout.build_h1_ledger(
        bins, cohort, roots, development_table, development_caches,
        root_file_sha256, development_table_sha256, route_context_sha256)
    binary = (result_dir / "frozen/wirehair_v2_bench").resolve(strict=True)
    cpus = tuple(contract.get("cpu_set", ()))
    thermal_policy = campaign.validate_thermal_policy(contract)
    taskset = holdout.frozen_taskset_path(contract)
    phase_record = holdout._load_phase_completion(
        result_dir, ledger, binary, set(cpus),
        thermal_policy=thermal_policy, taskset_path=taskset)
    phase_sha256 = verify_stem_hash_sidecar(
        result_dir / "h1/phase_complete.json")
    if phase_sha256 != sha256_bytes(
            campaign.canonical_json_bytes(phase_record)):
        die("H1 phase-completion semantic hash changed")

    index_sha256 = verify_named_hash_sidecar(
        result_dir / "h1/route_cache/index.json")
    controller_path = result_dir / "h1/controller_complete.json"
    controller = campaign.load_canonical_object(
        controller_path, "H1 controller completion")
    campaign.verify_sealed_record(
        controller,
        "wirehair.wh2.h12_preferred_attempt.h1_controller_complete.v1")
    verify_named_hash_sidecar(controller_path)
    expected_controller = h1_controller_complete_record(
        campaign,
        phase_completion_sha256=phase_sha256,
        analysis_complete_sha256=analysis_sha256,
        root_file_sha256=root_file_sha256,
        parent_table_sha256=development_table_sha256,
        frozen_table_sha256=h1_table_sha256,
        route_cache_index_sha256=index_sha256,
        route_context_sha256=route_context_sha256,
        input_routed_K=sum(
            value is not None for value in development_table.values()),
        retained_routed_K=sum(value is not None for value in h1_table.values()),
    )
    if controller != expected_controller:
        die("H1 controller-complete binding changed")

    expected_h1_files = {
        result_dir / "h1/ledgers/paired.json",
        result_dir / "h1/ledgers/paired.json.sha256",
        result_dir / "h1/phase_complete.json",
        result_dir / "h1/phase_complete.sha256",
        result_dir / "h1/analysis_complete.json",
        result_dir / "h1/analysis_complete.sha256",
        result_dir / "h1/frozen_table.json",
        result_dir / "h1/route_cache/index.json",
        result_dir / "h1/route_cache/index.json.sha256",
        controller_path,
        controller_path.with_suffix(".json.sha256"),
    }
    for binding in h1_caches.values():
        path = result_dir / binding.relative_path
        expected_h1_files.update((path, path.with_suffix(".csv.sha256")))
    for job in ledger.jobs:
        expected_h1_files.update(holdout.job_paths(result_dir, job).values())
    expected_h1_files = {
        path.resolve(strict=True) for path in expected_h1_files}
    actual_h1_files = exact_tree_files(result_dir, ("h1",))
    if actual_h1_files != expected_h1_files:
        die(
            "H1 terminal artifact inventory changed "
            f"(missing={len(expected_h1_files - actual_h1_files)}, "
            f"extra={len(actual_h1_files - expected_h1_files)})")


def verify_phase_ready_for_seal(
    result_dir: Path,
    spec: HoldoutPhaseSpec,
    prepare_record: dict[str, Any],
    contract: dict[str, Any],
) -> None:
    if spec.phase == "h1":
        verify_development_seal_readiness(
            result_dir, prepare_record, contract)
    elif spec.phase == "h2":
        verify_h1_seal_readiness(result_dir, contract)
    else:
        die("unknown holdout phase readiness contract")


def rooted_holdout_record(
    result_dir: Path,
    phase: str,
    contract: dict[str, Any],
) -> dict[str, Any]:
    spec = HOLDOUT_PHASES[phase]
    state = load_fixed_json(
        seal_state_path(result_dir, phase), f"{phase} controller state")
    if state.get("phase") != phase or state.get("status") != "ROOTED":
        die(f"{phase} holdout root has not been acquired")
    attempt, manifest, seal_record, seal_sha256 = load_seal_attempt(
        result_dir, spec, state)
    validate_loaded_seal(spec, manifest, seal_record, contract)
    validate_state_binding(state, spec, manifest, seal_record, seal_sha256)
    verify_manifest_bytes(result_dir, spec, manifest)
    verify_seal_publication(
        result_dir, state, attempt, seal_record, seal_sha256, contract)
    record = verify_rooted_record(
        result_dir, spec, contract, seal_sha256, seal_record)
    validate_rooted_state_binding(result_dir, spec, state)
    return record


def run_h1(args: argparse.Namespace) -> int:
    result_dir = args.result_dir.resolve(strict=True)
    _prepare, contract = verify_frozen_controller_runtime(result_dir)
    campaign = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_campaign",
        "campaign_module_sha256")
    holdout = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_holdout",
        "holdout_module_sha256")
    cohort, bins = load_frozen_cohort_bins(result_dir, campaign)
    route_context_sha256 = str(contract.get("route_context_sha256", ""))
    table, route_caches, parent_table_sha256 = \
        load_development_table_and_caches(
            result_dir, campaign, holdout, cohort, bins,
            route_context_sha256)
    rooted = rooted_holdout_record(result_dir, "h1", contract)
    roots = rooted.get("roots")
    if (not isinstance(roots, list) or len(roots) != H1_ROOT_COUNT or
            any(not isinstance(root, str) or
                not re.fullmatch(r"0x[0-9a-f]{16}", root) for root in roots)):
        die("H1 rooted seed ledger is malformed")
    root_file_sha256 = common.sha256_file(
        result_dir / HOLDOUT_PHASES["h1"].root_file)
    ledger = holdout.build_h1_ledger(
        bins, cohort, roots, table, route_caches,
        root_file_sha256, parent_table_sha256, route_context_sha256)
    binary = (result_dir / "frozen/wirehair_v2_bench").resolve(strict=True)
    cpus = tuple(contract.get("cpu_set", ()))
    runner = holdout.HoldoutRunner(
        result_dir, binary, cpus, contract.get("workers"),
        contract.get("timeout_seconds"))
    thermal = Path(str(contract.get("thermal", ""))).resolve(strict=True)
    runner.run_ledger(ledger, thermal)
    analysis = holdout.analyze_h1(
        result_dir, binary, ledger, cohort, table, set(cpus))
    frozen_table_record = analysis.get("frozen_table")
    if not isinstance(frozen_table_record, dict):
        die("H1 analysis omitted its frozen table")
    h1_table = holdout.parse_h1_table(frozen_table_record, strict_domain=True)

    h1_cache_rows: list[dict[str, Any]] = []
    for bin_index in range(BIN_COUNT):
        Ks = tuple(bins[bin_index])
        source_bytes = common.stable_bytes(
            result_dir / route_caches[bin_index].relative_path)
        source_records = holdout.parse_route_manifest(
            source_bytes, Ks, WIDTHS, {K: table[K] for K in Ks},
            route_context_sha256, True)
        relative = f"h1/route_cache/bin-{bin_index:03d}.csv"
        binding = holdout.publish_derived_route_cache(
            result_dir, relative, source_records, Ks, WIDTHS,
            {K: h1_table[K] for K in Ks}, route_context_sha256)
        h1_cache_rows.append({
            "bin": bin_index, "K_count": len(Ks),
            "path": binding.relative_path, "sha256": binding.sha256,
        })
    h1_table_path = result_dir / "h1/frozen_table.json"
    h1_table_sha256 = common.sha256_file(h1_table_path)
    index_record = campaign.sealed_record(
        "wirehair.wh2.h12_preferred_attempt.h1_route_cache_index.v1",
        {
            "frozen_table_sha256": h1_table_sha256,
            "route_context_sha256": route_context_sha256,
            "cache_count": len(h1_cache_rows), "caches": h1_cache_rows,
        },
    )
    index_path = result_dir / "h1/route_cache/index.json"
    index_sha256 = campaign.write_hashed_artifact(
        index_path, campaign.canonical_json_bytes(index_record))
    (validated_h1_table, _validated_h1_caches,
     validated_h1_table_sha256, analysis_sha256) = \
        load_h1_table_and_caches(
            result_dir, campaign, holdout, cohort, bins, table,
            parent_table_sha256, route_context_sha256)
    if (validated_h1_table != h1_table or
            validated_h1_table_sha256 != h1_table_sha256 or
            verify_named_hash_sidecar(index_path) != index_sha256):
        die("H1 terminal table/cache validation changed its published identity")
    phase_sha256 = verify_stem_hash_sidecar(
        result_dir / "h1/phase_complete.json")
    input_routed_K = sum(value is not None for value in table.values())
    retained_routed_K = sum(
        value is not None for value in h1_table.values())
    controller_record = h1_controller_complete_record(
        campaign,
        phase_completion_sha256=phase_sha256,
        analysis_complete_sha256=analysis_sha256,
        root_file_sha256=root_file_sha256,
        parent_table_sha256=parent_table_sha256,
        frozen_table_sha256=h1_table_sha256,
        route_cache_index_sha256=index_sha256,
        route_context_sha256=route_context_sha256,
        input_routed_K=input_routed_K,
        retained_routed_K=retained_routed_K,
    )
    controller_path = result_dir / "h1/controller_complete.json"
    controller_sha256 = campaign.write_hashed_artifact(
        controller_path, campaign.canonical_json_bytes(controller_record))
    complete = {
        "schema": "wirehair.wh2.h12_preferred_attempt.h1_complete.v1",
        "accepted": analysis.get("accepted"),
        "input_routed_K": input_routed_K,
        "retained_routed_K": retained_routed_K,
        "frozen_table_sha256": h1_table_sha256,
        "route_cache_index_sha256": index_sha256,
        "controller_complete_sha256": controller_sha256,
    }
    print(canonical_json(complete))
    return 0


def run_h2(args: argparse.Namespace) -> int:
    result_dir = args.result_dir.resolve(strict=True)
    _prepare, contract = verify_frozen_controller_runtime(result_dir)
    campaign = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_campaign",
        "campaign_module_sha256")
    holdout = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_holdout",
        "holdout_module_sha256")
    all_K, all_bins = load_frozen_allk_bins(result_dir, campaign)
    route_context_sha256 = str(contract.get("route_context_sha256", ""))
    h1_table_path = result_dir / "h1/frozen_table.json"
    h1_record = campaign.load_canonical_object(
        h1_table_path, "H1 frozen table")
    active_table = holdout.parse_h1_table(h1_record, strict_domain=True)
    table = holdout.expand_h1_table_all_k(h1_record)
    if tuple(table) != all_K or any(
            table[K] is not None for K in all_K if K not in active_table):
        die("H1 all-K expansion changed non-active controls")
    parent_table_sha256 = common.sha256_file(h1_table_path)
    rooted = rooted_holdout_record(result_dir, "h2", contract)
    roots = rooted.get("roots")
    if (not isinstance(roots, list) or len(roots) != H2_ROOT_COUNT or
            any(not isinstance(root, str) or
                not re.fullmatch(r"0x[0-9a-f]{16}", root) for root in roots)):
        die("H2 rooted seed ledger is malformed")
    root_file_sha256 = common.sha256_file(
        result_dir / HOLDOUT_PHASES["h2"].root_file)
    route_ledger = holdout.build_h2_route_ledger(
        all_bins, roots, table, root_file_sha256,
        parent_table_sha256, route_context_sha256)
    binary = (result_dir / "frozen/wirehair_v2_bench").resolve(strict=True)
    cpus = tuple(contract.get("cpu_set", ()))
    runner = holdout.HoldoutRunner(
        result_dir, binary, cpus, contract.get("workers"),
        contract.get("timeout_seconds"))
    thermal = Path(str(contract.get("thermal", ""))).resolve(strict=True)
    runner.run_ledger(route_ledger, thermal)
    route_caches = holdout.derive_h2_bb64_route_caches(
        result_dir, binary, route_ledger, table, set(cpus))
    paired_ledger = holdout.build_h2_paired_ledger(
        all_bins, roots, table, route_caches, root_file_sha256,
        parent_table_sha256, route_context_sha256)
    runner.run_ledger(paired_ledger, thermal)
    analysis = holdout.analyze_h2(
        result_dir, binary, route_ledger, paired_ledger, table, set(cpus))
    complete = {
        "schema": "wirehair.wh2.h12_preferred_attempt.h2_complete.v1",
        "accepted": analysis.get("accepted"),
        "routed_K": sum(value is not None for value in table.values()),
        "control_K": sum(value is None for value in table.values()),
        "phase_complete_sha256": common.sha256_file(
            result_dir / "h2/phase_complete.json"),
        "analysis_complete_sha256": common.sha256_file(
            result_dir / "h2/analysis_complete.json"),
    }
    print(canonical_json(complete))
    return 0


def timing_panel_specs(
    timing: Any,
    sample: Sequence[int],
    table: dict[int, int | None],
) -> tuple[Any, ...]:
    if (tuple(timing.WIDTHS) != (64, 1280, 4096) or
            tuple(timing.SCHEDULES) != SCHEDULES):
        die("PERMANENT: timing module grid differs from the frozen protocol")
    specs: list[Any] = []
    for K in sample:
        preferred = table.get(K)
        if not isinstance(preferred, int) or isinstance(preferred, bool):
            die("timing sample contains an unrouted K")
        for width in timing.WIDTHS:
            for schedule in timing.SCHEDULES:
                controller_seed = timing_panel_seed(
                    "solve", K, width, schedule)
                module_seed = timing.timing_panel_seed(
                    "solve", K, width, schedule)
                if controller_seed != module_seed:
                    die("PERMANENT: timing seed derivations disagree")
                specs.append(timing.TimingPanelSpec(
                    K, width, preferred, "solve", controller_seed, schedule))
        for width in timing.WIDTHS:
            controller_seed = timing_panel_seed(
                "setup", K, width, TIMING_SETUP_SCHEDULE)
            module_seed = timing.timing_panel_seed(
                "setup", K, width, TIMING_SETUP_SCHEDULE)
            if controller_seed != module_seed:
                die("PERMANENT: timing setup-seed derivations disagree")
            specs.append(timing.TimingPanelSpec(
                K, width, preferred, "setup", controller_seed,
                TIMING_SETUP_SCHEDULE))
    if len(specs) != 12 * len(sample) or len({spec.key() for spec in specs}) != len(specs):
        die("timing panel specification coverage changed")
    return tuple(specs)


def timing_panel_paths(
    timing_dir: Path,
    spec: Any,
) -> tuple[Path, Path]:
    stem = f"K-{spec.K:05d}-bb-{spec.block_bytes}-{spec.schedule}"
    return (
        timing_dir / "panels" / spec.metric / stem,
        timing_dir / "evidence" / spec.metric / stem,
    )


def verify_timing_evidence_files(
    timing_dir: Path,
    evidence_dir: Path,
    result: Any,
) -> set[Path]:
    if evidence_dir.is_symlink() or not evidence_dir.is_dir():
        die("timing panel lacks its raw evidence directory")
    expected: dict[str, str] = {}
    for attempt in result.attempts:
        suffix = (
            "full" if attempt.cycle_index is None else
            f"cycle-{attempt.cycle_index}")
        stem = f"attempt-{attempt.attempt_index:02d}-{suffix}"
        expected[stem + ".thermal.csv"] = \
            attempt.environment.thermal_interval_sha256
        expected[stem + ".performance.bin"] = \
            attempt.environment.performance_interval_sha256
    actual = {path.name: path for path in evidence_dir.iterdir()}
    if set(actual) != set(expected):
        die("timing raw-evidence inventory is not exact")
    verified: set[Path] = set()
    for name, digest in expected.items():
        path = actual[name]
        if (path.is_symlink() or not path.is_file() or
                common.sha256_file(path) != digest):
            die("timing raw evidence changed after its attempt was sealed")
        try:
            path.relative_to(timing_dir)
        except ValueError:
            die("timing raw evidence escaped the timing directory")
        verified.add(path.resolve(strict=True))
    return verified


def publish_combined_timing_cache(
    result_dir: Path,
    timing_dir: Path,
    campaign: Any,
    holdout: Any,
    sample: Sequence[int],
    table: dict[int, int | None],
    bins: dict[int, tuple[int, ...]],
    bindings: dict[int, Any],
    route_context_sha256: str,
) -> tuple[Path, str]:
    selected = set(sample)
    source_records: list[Any] = []
    for bin_index in range(BIN_COUNT):
        Ks = bins[bin_index]
        binding = bindings[bin_index]
        records = holdout.parse_route_manifest(
            common.stable_bytes(result_dir / binding.relative_path),
            Ks, WIDTHS, {K: table[K] for K in Ks},
            route_context_sha256, True)
        source_records.extend(record for record in records if record.K in selected)
    cache_bytes = holdout.canonical_selected_route_manifest(
        source_records, tuple(sample), WIDTHS,
        {K: table[K] for K in sample}, route_context_sha256)
    cache_path = timing_dir / "route_cache.csv"
    cache_sha256 = campaign.write_hashed_artifact(cache_path, cache_bytes)
    holdout.parse_route_manifest(
        common.stable_bytes(cache_path), tuple(sample), WIDTHS,
        {K: table[K] for K in sample}, route_context_sha256, True)
    if verify_named_hash_sidecar(cache_path) != cache_sha256:
        die("combined timing route-cache binding changed")
    return cache_path.resolve(strict=True), cache_sha256


def run_timing(args: argparse.Namespace) -> int:
    """Run or strictly resume the final isolated H1-table timing gate."""
    result_dir = args.result_dir.resolve(strict=True)
    _prepare, contract = verify_frozen_controller_runtime(result_dir)
    campaign = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_campaign",
        "campaign_module_sha256")
    holdout = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_holdout",
        "holdout_module_sha256")
    timing = load_frozen_python_module(
        result_dir, "wh2_preferred_attempt_timing",
        "timing_module_sha256")
    cohort, bins = load_frozen_cohort_bins(result_dir, campaign)
    route_context_sha256 = str(contract.get("route_context_sha256", ""))
    development_table, _development_caches, development_table_sha256 = \
        load_development_table_and_caches(
            result_dir, campaign, holdout, cohort, bins,
            route_context_sha256)
    h1_table, h1_caches, h1_table_sha256, h1_analysis_sha256 = \
        load_h1_table_and_caches(
            result_dir, campaign, holdout, cohort, bins,
            development_table, development_table_sha256,
            route_context_sha256)

    # Commit and reveal the independent H2 domain before any timing result can
    # influence whether that final recovery holdout is acquired.
    h2_root = rooted_holdout_record(result_dir, "h2", contract)
    h2_root_path = result_dir / HOLDOUT_PHASES["h2"].root_file
    h2_root_sha256 = common.sha256_file(h2_root_path.resolve(strict=True))
    routed = tuple(K for K in cohort if h1_table[K] is not None)
    timing_dir = result_dir / "timing"
    timing_dir.mkdir(mode=0o700, exist_ok=True)
    if timing_dir.is_symlink() or not timing_dir.is_dir():
        die("timing result directory is not a regular directory")
    if len(routed) < 5:
        rejection = campaign.sealed_record(
            "wirehair.wh2.h12_preferred_attempt.timing_rejection.v1",
            {
                "accepted": False,
                "reason": "fewer_than_five_routed_K",
                "routed_K": len(routed),
                "h1_table_sha256": h1_table_sha256,
                "h1_analysis_sha256": h1_analysis_sha256,
                "h2_root_file_sha256": h2_root_sha256,
                "h2_seal_record_sha256": h2_root["seal_record_sha256"],
            },
        )
        path = timing_dir / "rejection.json"
        digest = campaign.write_hashed_artifact(
            path, campaign.canonical_json_bytes(rejection))
        expected = {
            path.resolve(strict=True),
            path.with_suffix(".json.sha256").resolve(strict=True),
        }
        if exact_tree_files(result_dir, ("timing",)) != expected:
            die("insufficient-sample timing rejection tree is not exact")
        print(canonical_json({
            "schema": rejection["schema"], "accepted": False,
            "report_sha256": digest, "routed_K": len(routed),
        }))
        return 0

    try:
        sample = timing.select_timing_sample(routed)
    except Exception as error:
        die(f"timing sample selection failed: {error}")
    specs = timing_panel_specs(timing, sample, h1_table)
    cache_path, cache_sha256 = publish_combined_timing_cache(
        result_dir, timing_dir, campaign, holdout, sample, h1_table,
        bins, h1_caches, route_context_sha256)
    sample_record = campaign.sealed_record(
        "wirehair.wh2.h12_preferred_attempt.timing_sample.v1",
        {
            "protocol_sha256": contract["protocol_sha256"],
            "route_context_sha256": route_context_sha256,
            "h1_table_sha256": h1_table_sha256,
            "h1_analysis_sha256": h1_analysis_sha256,
            "h2_root_file_sha256": h2_root_sha256,
            "h2_seal_record_sha256": h2_root["seal_record_sha256"],
            "execution_core": TIMING_CORE,
            "setup_schedule": TIMING_SETUP_SCHEDULE,
            "sample": list(sample), "S": len(sample),
            "route_cache_sha256": cache_sha256,
            "panel_count": len(specs),
            "panels": [
                {
                    "metric": spec.metric, "K": spec.K,
                    "block_bytes": spec.block_bytes,
                    "schedule": spec.schedule,
                    "preferred_attempt": spec.preferred_attempt,
                    "seed": spec.seed,
                }
                for spec in specs
            ],
        },
    )
    sample_path = timing_dir / "sample.json"
    sample_sha256 = campaign.write_hashed_artifact(
        sample_path, campaign.canonical_json_bytes(sample_record))

    binary = (result_dir / "frozen/wirehair_v2_bench").resolve(strict=True)
    thermal_log = Path(str(contract.get("thermal", ""))).resolve(strict=True)
    try:
        host = TimingHostSession(contract, timing, TIMING_CORE)
        config = timing.TimingRunnerConfig(
            binary=binary,
            binary_sha256=str(contract.get("binary_sha256", "")),
            route_cache=cache_path,
            route_cache_sha256=cache_sha256,
            route_context_sha256=route_context_sha256,
            isolation=host.isolation,
            launcher=host.launcher_spec(),
            timeout_seconds=float(contract.get("timeout_seconds", 0)),
        )
    except Exception as error:
        die(f"timing host/configuration validation failed: {error}")

    panels: list[Any] = []
    execution_files: set[Path] = {
        cache_path,
        cache_path.with_suffix(".csv.sha256").resolve(strict=True),
        sample_path.resolve(strict=True),
        sample_path.with_suffix(".json.sha256").resolve(strict=True),
    }
    missing_panel = any(
        not timing_panel_paths(timing_dir, spec)[0].exists()
        for spec in specs)

    def load_or_run_panels() -> None:
        for spec in specs:
            panel_dir, evidence_dir = timing_panel_paths(timing_dir, spec)
            probe = LinuxTimingEvidenceProbe(
                timing, host, thermal_log, evidence_dir)
            try:
                result = timing.run_or_resume_timing_panel(
                    spec, config, probe, panel_dir)
            except Exception as error:
                die(
                    "timing panel failed for "
                    f"{spec.metric}/K={spec.K}/bb={spec.block_bytes}/"
                    f"{spec.schedule}: {error}")
            panels.append(result)
            for path in panel_dir.iterdir():
                if path.is_symlink() or not path.is_file():
                    die("timing panel contains a nonregular artifact")
                execution_files.add(path.resolve(strict=True))
            execution_files.update(verify_timing_evidence_files(
                timing_dir, evidence_dir, result))

    if missing_panel:
        with host:
            if host.launcher() != config.launcher:
                die("timing launcher changed when isolation became live")
            load_or_run_panels()
    else:
        load_or_run_panels()

    phase_sha256 = write_exact_manifest_once(
        result_dir, timing_dir / "phase_complete.sha256",
        execution_files, campaign)
    try:
        analysis = timing.analyze_timing(
            sample, panels, setup_schedule=TIMING_SETUP_SCHEDULE)
    except Exception as error:
        die(f"timing exact analysis failed: {error}")
    report = campaign.sealed_record(
        "wirehair.wh2.h12_preferred_attempt.timing_complete.v1",
        {
            "phase_seal_sha256": phase_sha256,
            "sample_sha256": sample_sha256,
            "route_cache_sha256": cache_sha256,
            "h1_table_sha256": h1_table_sha256,
            "h1_analysis_sha256": h1_analysis_sha256,
            "h2_root_file_sha256": h2_root_sha256,
            "analysis": analysis,
            "accepted": analysis.get("accepted"),
        },
    )
    report_path = timing_dir / "report.json"
    report_sha256 = campaign.write_hashed_artifact(
        report_path, campaign.canonical_json_bytes(report))
    analysis_files = set(execution_files)
    analysis_files.update((
        (timing_dir / "phase_complete.sha256").resolve(strict=True),
        report_path.resolve(strict=True),
        report_path.with_suffix(".json.sha256").resolve(strict=True),
    ))
    analysis_sha256 = write_exact_manifest_once(
        result_dir, timing_dir / "analysis_complete.sha256",
        analysis_files, campaign)
    expected_tree = set(analysis_files)
    expected_tree.add(
        (timing_dir / "analysis_complete.sha256").resolve(strict=True))
    if exact_tree_files(result_dir, ("timing",)) != expected_tree:
        die("completed timing result tree has unexpected artifacts")
    complete = {
        "schema": "wirehair.wh2.h12_preferred_attempt.timing_run.v1",
        "accepted": analysis.get("accepted"), "S": len(sample),
        "phase_seal_sha256": phase_sha256,
        "analysis_seal_sha256": analysis_sha256,
        "report_sha256": report_sha256,
    }
    print(canonical_json(complete))
    return 0


def plan(args: argparse.Namespace) -> int:
    protocol, ledgers = derive_plan(args.source_result)
    report = {
        "schema": protocol["schema"],
        "source_commit": SOURCE_COMMIT,
        "protocol_sha256": sha256_bytes(common.json_bytes(protocol)),
        "ledgers": {
            name: {"bytes": len(data), "sha256": sha256_bytes(data)}
            for name, data in ledgers.items()
        },
        "job_totals": protocol["job_totals"],
        "development": {
            key: protocol["development"][key]
            for key in ("candidate_cells", "full_control_cube_cells",
                        "logical_cells", "candidate_jobs", "control_jobs", "jobs")
        },
        "h1": {
            "logical_cells": protocol["h1"]["logical_cells"],
            "jobs": protocol["h1"]["jobs"],
        },
        "h2": {
            "logical_cells": protocol["h2"]["logical_cells"],
            "jobs": protocol["h2"]["jobs"],
        },
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    plan_parser = subparsers.add_parser("plan", help="derive and audit the frozen plan")
    plan_parser.add_argument("--source-result", type=Path, required=True)
    prepare_parser = subparsers.add_parser(
        "prepare", help="freeze immutable inputs without launching recovery work")
    prepare_parser.add_argument("--source-result", type=Path, required=True)
    prepare_parser.add_argument("--binary", type=Path, required=True)
    prepare_parser.add_argument("--thermal", type=Path, required=True)
    prepare_parser.add_argument("--result-dir", type=Path, required=True)
    prepare_parser.add_argument("--workers", type=int, default=128)
    prepare_parser.add_argument("--build-workers", type=int, default=128)
    prepare_parser.add_argument("--timeout", type=float, default=1800.0)
    prepare_parser.add_argument("--acknowledge", required=True)
    development_parser = subparsers.add_parser(
        "run-development",
        help="resume control and R1..R4 from the exact frozen controller",
    )
    development_parser.add_argument("--result-dir", type=Path, required=True)
    h1_parser = subparsers.add_parser(
        "run-h1", help="resume the exact rooted H1 pruning campaign")
    h1_parser.add_argument("--result-dir", type=Path, required=True)
    h2_parser = subparsers.add_parser(
        "run-h2", help="resume the exact rooted all-K H2 campaign")
    h2_parser.add_argument("--result-dir", type=Path, required=True)
    timing_parser = subparsers.add_parser(
        "run-timing", help="run the exact isolated H1-table timing gate")
    timing_parser.add_argument("--result-dir", type=Path, required=True)
    seal_parser = subparsers.add_parser(
        "seal-holdout",
        help="seal a completed table and publish its future-round GitHub tag",
    )
    seal_parser.add_argument("--result-dir", type=Path, required=True)
    seal_parser.add_argument("--phase", choices=tuple(HOLDOUT_PHASES), required=True)
    acquire_parser = subparsers.add_parser(
        "acquire-holdout",
        help="resume the exact sealed round until a verified holdout root exists",
    )
    acquire_parser.add_argument("--result-dir", type=Path, required=True)
    acquire_parser.add_argument(
        "--phase", choices=tuple(HOLDOUT_PHASES), required=True)
    acquire_parser.add_argument("--max-wait-seconds", type=float, default=900.0)
    acquire_parser.add_argument("--no-wait", action="store_true")
    args = parser.parse_args(argv)
    if hasattr(args, "timeout") and (
            not math.isfinite(args.timeout) or args.timeout <= 0):
        parser.error("--timeout must be finite and positive")
    if hasattr(args, "max_wait_seconds") and (
            not math.isfinite(args.max_wait_seconds) or
            args.max_wait_seconds <= 0 or args.max_wait_seconds > 900):
        parser.error("--max-wait-seconds must be in (0,900]")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.command == "plan":
        return plan(args)
    if args.command == "prepare":
        return prepare(args)
    if args.command == "run-development":
        return run_development(args)
    if args.command == "run-h1":
        return run_h1(args)
    if args.command == "run-h2":
        return run_h2(args)
    if args.command == "run-timing":
        return run_timing(args)
    if args.command == "seal-holdout":
        return seal_holdout(args)
    if args.command == "acquire-holdout":
        return acquire_holdout(args)
    die(f"unknown command {args.command}")
    return 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CampaignError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
