#!/usr/bin/env python3
"""Shared, fail-closed access to the native WH2 peel benchmark.

The native codec is the authority for the shipped peel distribution.  Do not
reconstruct it from an analytic soliton formula here: the production generator
has a fixed degree-one prefix through K=2048 and no degree-one prefix above it,
and its fixed-point thresholds are not the same as ``1 / K``.
"""

import hashlib
import json
import math
import os
import platform
import stat
import subprocess
import sys
import tempfile
from dataclasses import dataclass


PEEL_TABLE_SCHEMA = "wirehair-v2-peel-table/v3"
NATIVE_PMF_PROTOCOL = "wirehair-v2-bench:peelpmf:dispatch-v1:v2"
NATIVE_COMPARE_PROTOCOL = "wirehair-v2-bench:compare:wh2-target:v3"
SEED_DERIVATION_PROTOCOL = "blake2b64:wh2-peel-v2:split-construction-loss"
COMPLETION_REGIME_PROTOCOL = "wirehair-v2-dispatch-v1-completion/v2"
TARGET_PROFILE = "dispatch-v1"
TARGET_SEED_POLICY = "raw"
TARGET_RECEIPT_SEED_POLICY = "raw-uniform-v1"
TARGET_CONTRACT = {
    "target_profile": TARGET_PROFILE,
    "contract_id": "a98c37c23ee7feae",
    "contract_sha256":
        "a98c37c23ee7feae4171ff10627f660f705db6b7aae9268f85617ce86396583c",
    "precode_contract": 5,
    "packet_contract": 4,
    "architecture": "smallband100-d4",
    "completion": "mixed",
    "dense_rows": 4,
    "heavy_rows": 12,
    "gf256_rows": 10,
    "gf16_rows": 2,
    "period": 244,
    "geometry": "frozen",
    "residue_schedule": "constant",
    "residue_skew": 0,
    "mix_count": 3,
    "pmf_encoding": "wirehair-v2-peel-spec-v1",
}
TARGET_MEASUREMENT_POLICY = {
    "seed_policy": TARGET_RECEIPT_SEED_POLICY,
    "attempt_policy": "single",
    "seed_attempt": "0",
    "seed_tables": "none",
    "fixups": "none",
    "grouped_gf256_rows": "0",
    "independent_extension": "0",
    "heavy_family": "periodic",
    "dense_identity_corner": "0",
    "dense_two_anchor": "0",
    "packet_seed_multiplier": "1",
    "packet_seed_avalanche": "0",
    "loss": "iid-drop",
    "loss_trace": "common-id-v2",
}
# The current native receipt says ``loss=iid-drop``.  Until that field is
# versioned into a schedule-specific loss model, accepting a non-IID schedule
# here would publish a false receipt.
TARGET_SCHEDULES = {"iid"}
STOCK_PMF_DIGEST = hashlib.sha256(b"stock").hexdigest()
PROXY_COST_MODEL = "embedded-es-cost-model/raw-calibration-unavailable"
PROXY_MEASURE_REGIME = {
    "solve_block_bytes": 0,
    "cost_model_block_bytes": 1280,
    "cost_model_verified": 0,
    "band_tracking_x": 1,
    "loss": "0.100000",
    "seed_base": 55,
    "completion": "mixed",
    "geometry": "frozen",
    "period": 244,
    "gf16_rows": 0,
}
PROXY_ORDERING_PROTOCOL = "fail-rate-then-pred-ns/v1"
SEARCH_BOX_PROTOCOL = "native-mean:[0,min(K,max(40,4*mean))]/v1"
PROXY_K_LADDER = (
    2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384,
    512, 768, 1024, 1536, 2048, 4096, 8192, 16384, 32768, 64000,
)
# These are the complete production mixed-completion settings reported by the
# native benchmark.  Search tools do not pass test-hook overrides for them:
# if a production default changes, the receipt must get a deliberate schema
# update rather than silently comparing a different codec under the old name.
PRODUCTION_COMPLETION_REGIME = {
    "precode": "1",
    "precode_cache": "0",
    "precode_profile": "mixed",
    "encoder_cache": "0",
    "decoder_cache": "0",
    "mixed_mix_count": "3",
    "mixed_period": "244",
    "mixed_gf256_rows": "10",
    "mixed_gf16_rows": "2",
    "mixed_geometry": "frozen",
    "mixed_residue_skew": "0",
    "mixed_residue_schedule": "constant",
    "mixed_residue_hash_seed": "0x0",
    "mixed_residue_hash_keyed": "0",
    "mixed_independent_extension_residues": "0",
    "mixed_grouped_gf256_rows": "0",
    "mixed_grouped_gf256_row_mask": "0x0",
    "mixed_grouped_gf256_hash_seed": "0x0",
    "mixed_grouped_final_h_a_columns": "0",
    "mixed_residue_buckets_requested": "auto",
    "packet_row_seed_multiplier": "0x1",
    "packet_row_seed_avalanche": "0",
    "precode_profile_handoff": "encoder-selected-v1",
}
# ValidatePrecodeParams() bounds K + S + D2 + H to a uint16 span.  The exact
# dispatch-v1 target fixes D2 at four rows and its completion height to the
# GF(256) plus GF(2^16) rows named in the receipt above.  Native PMF metadata
# comes from that contract, so a larger S could not have been produced by a
# valid target expansion.
_PRODUCTION_DENSE_ROWS = int(TARGET_CONTRACT["dense_rows"])
_PRODUCTION_TOTAL_SPAN_MAX = 0xffff
_PRODUCTION_MIXED_HEAVY_ROWS = (
    int(PRODUCTION_COMPLETION_REGIME["mixed_gf256_rows"]) +
    int(PRODUCTION_COMPLETION_REGIME["mixed_gf16_rows"])
)
# Exact WirehairTools::GetDenseCount() data used by dispatch-v1's
# SmallBandStaircaseCount().  Keep the legacy values here instead of trusting
# a self-reported native staircase: table documents are security boundaries,
# and a forged receipt must not be able to redefine the equation geometry.
_LEGACY_TINY_DENSE_COUNTS = (
    0, 0, 2, 3, 3, 5, 6, 6, 6, 7, 9, 10, 10, 10, 12, 14,
    13, 14, 12, 12, 15, 16, 21, 14, 14, 13, 18, 21, 22, 21, 13, 22,
    13, 24, 14, 17, 16, 24, 30, 26, 24, 18, 15, 15, 24, 18, 21, 17,
    14, 16, 21, 18, 17, 22, 25, 20, 17, 18, 21, 18, 23, 20, 19, 23,
    19,
)
_LEGACY_DENSE_POINTS = (
    (2048, 52),
    (2618, 54),
    (2826, 60),
    (3725, 62),
    (3962, 67),
    (4277, 65),
    (4547, 60),
    (5065, 64),
    (5224, 76),
    (5642, 66),
    (5909, 71),
    (6285, 76),
    (6583, 66),
    (6895, 72),
    (7448, 69),
    (7682, 76),
    (8046, 78),
    (8558, 76),
    (8963, 73),
    (9389, 81),
    (10143, 86),
    (11129, 94),
    (12593, 99),
    (12988, 105),
    (14032, 108),
    (14473, 114),
    (15397, 110),
    (16636, 113),
    (17698, 118),
    (18828, 123),
    (19420, 127),
    (20343, 136),
    (21979, 139),
    (23024, 150),
    (24119, 156),
    (25659, 162),
    (27298, 173),
    (29042, 176),
    (30898, 183),
    (31870, 190),
    (33906, 200),
    (35519, 211),
    (37208, 220),
    (38978, 234),
    (40205, 253),
    (42776, 297),
    (44122, 320),
    (45511, 336),
    (46944, 357),
    (48421, 373),
    (49177, 376),
    (50725, 380),
    (52321, 391),
    (53968, 388),
    (54811, 382),
    # This duplicate is intentional and mirrors the native graph exactly.
    (54811, 382),
    (55667, 372),
    (57419, 362),
    (58316, 356),
    (60152, 347),
    (61091, 337),
    (62045, 334),
    (63014, 340),
    (64000, 345),
)
# The native compare parser stores --trials in uint32_t and parses --bb-list
# through a positive int.  The production mixed arm additionally requires an
# even payload width for its GF(2^16) rows.
_COMPARE_TRIALS_MAX = 0xffffffff
_COMPARE_BLOCK_BYTES_MAX = 0x7fffffff
PRODUCTION_COMPARE_ARM = {
    "target_profile": TARGET_PROFILE,
    "seed_policy": TARGET_SEED_POLICY,
}
_COMPARE_BANNER_ARM = {
    "v2_profile": "base",
    "peel_candidates": "16",
    "peel_trials": "3",
    "auto_trials": "8",
    "auto_min_delta": "0.1000",
    "tune_seed": "0x9a7e11a",
    "auto_seed": "0xa570ca1",
    "dense_override": "0",
    "dense_delta": "0",
    "dense_candidate": "0",
}
RECOVERY_METRIC_SCOPE = {
    "fail": "all-trials",
    "overhead": "successful-trials-only",
    "throughput": "successful-trials-only",
}
RECOVERY_METRIC_FIELDS = {
    "construction_seed", "loss_seed", "target_receipt", "fail", "oh_mean",
    "OH_sd", "OH50", "OH95", "OH99", "OH_max", "decode_mbps",
}
COMPARE_COLUMNS = (
    "codec", "bb", "trials", "fail", "N_mean", "OH_mean", "OH_sd",
    "OH50", "OH95", "OH99", "OH_max", "create_MBps", "encode_MBps",
    "decode_MBps", "recover_MBps",
)
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SOURCE_EXTENSIONS = (".c", ".cc", ".cmake", ".cpp", ".h", ".hpp", ".inc", ".py")


class MeasurementError(RuntimeError):
    """The benchmark failed or emitted an unusable receipt."""


def valid_loss_rate(value):
    """Return whether value has the native benchmark's canonical loss form."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return False
    try:
        numeric = float(value)
    except (OverflowError, TypeError, ValueError):
        return False
    return (
        math.isfinite(numeric) and
        0.0 <= numeric <= 0.99 and
        not (numeric == 0.0 and math.copysign(1.0, numeric) < 0.0)
    )


def _canonical_staircase_scale_spec(value):
    """Return the native receipt spelling for a validated numeric scale."""
    numeric = float(value)
    if numeric == 0.0:
        numeric = 0.0
    return f"{numeric:.17g}"


def _target_mean_spec(
        block_count, staircase, source_hits, staircase_scale):
    """Return native's exact printed mean staircase-row degree."""
    if staircase_scale == "unset":
        edge_count = block_count * min(source_hits, staircase)
    else:
        scale = float(staircase_scale)
        # MakeStaircaseDegreeMixture() uses llround(scale * S).  Its input
        # domain is nonnegative, where llround is floor(x + 0.5).
        edge_count = math.floor(scale * staircase + 0.5)
        edge_count = min(edge_count, block_count * staircase)
    return f"{edge_count / staircase:.17g}"


@dataclass(frozen=True)
class StockProfile:
    """Exact dispatch-v1 peel distribution and construction contract."""

    block_count: int
    target_profile: str
    contract_id: str
    contract_sha256: str
    precode_contract: int
    packet_contract: int
    architecture: str
    staircase: int
    dense_rows: int
    heavy_rows: int
    source_hits: int
    completion: str
    gf256_rows: int
    gf16_rows: int
    period: int
    geometry: str
    residue_schedule: str
    residue_skew: int
    mix_count: int
    target_mean: float
    native_pmf_sha256: str
    pmf_encoding: str
    pmf: tuple

    def as_dict(self):
        return {
            "block_count": self.block_count,
            "target_profile": self.target_profile,
            "contract_id": self.contract_id,
            "contract_sha256": self.contract_sha256,
            "precode_contract": self.precode_contract,
            "packet_contract": self.packet_contract,
            "architecture": self.architecture,
            "staircase": self.staircase,
            "dense_rows": self.dense_rows,
            "heavy_rows": self.heavy_rows,
            "source_hits": self.source_hits,
            "completion": self.completion,
            "gf256_rows": self.gf256_rows,
            "gf16_rows": self.gf16_rows,
            "period": self.period,
            "geometry": self.geometry,
            "residue_schedule": self.residue_schedule,
            "residue_skew": self.residue_skew,
            "mix_count": self.mix_count,
            "target_mean": self.target_mean,
            "native_pmf_sha256": self.native_pmf_sha256,
            "pmf_encoding": self.pmf_encoding,
            "pmf_sha256": pmf_sha256(self.pmf),
            "pmf": list(self.pmf),
        }


@dataclass(frozen=True)
class RecoveryMetrics:
    """One exact raw-seed dispatch-v1 recovery measurement."""

    construction_seed: int
    loss_seed: int
    target_receipt: dict
    fail: int
    oh_mean: float
    oh_sd: float
    oh50: float
    oh95: float
    oh99: float
    oh_max: float
    decode_mbps: float

    def goodput(self, block_count):
        if self.fail != 0:
            return 0.0
        return self.decode_mbps * block_count / (block_count + self.oh_mean)

    def as_dict(self):
        return {
            "construction_seed": self.construction_seed,
            "loss_seed": self.loss_seed,
            "target_receipt": dict(self.target_receipt),
            "fail": self.fail,
            "oh_mean": self.oh_mean,
            "OH_sd": self.oh_sd,
            "OH50": self.oh50,
            "OH95": self.oh95,
            "OH99": self.oh99,
            "OH_max": self.oh_max,
            "decode_mbps": self.decode_mbps,
        }


def make_search_receipt(
        metrics, *, mode, goodput, trials, block_bytes, search_kind,
        construction_seed_base, loss_seed_base, seed_domain, coordinates,
        peel_pmf, shipped_control, shipped_goodput, context=None):
    """Create a complete selected-arm receipt for a search result."""
    receipt = metrics.as_dict()
    receipt.update({
        "mode": mode,
        "goodput": goodput,
        "trials": trials,
        "block_bytes": block_bytes,
        "search_kind": search_kind,
        "construction_seed_base": construction_seed_base,
        "loss_seed_base": loss_seed_base,
        "seed_domain": seed_domain,
        "coordinates": dict(coordinates),
        "peel_pmf": list(peel_pmf),
        "peel_pmf_sha256": pmf_sha256(peel_pmf),
        "shipped_control": shipped_control.as_dict(),
        "shipped_goodput": shipped_goodput,
        "context": context or {},
    })
    return receipt


def isolated_codec_env(overrides=None):
    """Return an environment with every inherited WH2 test hook removed."""
    env = {
        key: value for key, value in os.environ.items()
        if not key.startswith("WIREHAIR_V2_")
    }
    if overrides:
        env.update(overrides)
    return env


def _run_checked(command, env):
    try:
        result = subprocess.run(
            command, capture_output=True, text=True, env=env)
    except OSError as error:
        raise MeasurementError(
            f"could not execute benchmark {command[0]!r}: {error}")
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or "no output"
        raise MeasurementError(
            f"benchmark exited {result.returncode}: {' '.join(command)}: "
            f"{detail}")
    return result.stdout


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def file_sha256(path):
    """Return a file content identity as a lowercase SHA-256 digest."""
    try:
        return _sha256_file(path)
    except OSError as error:
        raise MeasurementError(f"could not hash {path!r}: {error}")


def _benchmark_identity(bench):
    """Hash a stable regular-file snapshot, detecting concurrent replacement."""
    path = os.path.realpath(bench)
    try:
        with open(path, "rb") as executable:
            before = os.fstat(executable.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise OSError("benchmark is not a regular file")
            digest = hashlib.sha256()
            while True:
                chunk = executable.read(1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
            after = os.fstat(executable.fileno())
        if not stat.S_ISREG(after.st_mode):
            raise OSError("benchmark is not a regular file")
    except OSError as error:
        raise MeasurementError(
            f"could not identify benchmark {bench!r}: {error}")
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns")
    if any(getattr(before, name) != getattr(after, name)
           for name in stable_fields):
        raise MeasurementError(
            f"benchmark changed while it was being identified: {path!r}")
    return {
        "path": path,
        "sha256": digest.hexdigest(),
        "size": before.st_size,
    }


def _artifact_identity(bench, generator):
    """Return immutable executable and source-state identities."""
    benchmark = _benchmark_identity(bench)
    generator_path = (
        generator if os.path.isabs(generator)
        else os.path.join(REPO_ROOT, generator)
    )
    try:
        generator_sha256 = _sha256_file(generator_path)
    except OSError as error:
        raise MeasurementError(
            f"could not identify generator source: {error}")

    try:
        git = subprocess.run(
            ["git", "-C", REPO_ROOT, "rev-parse", "HEAD"],
            capture_output=True, text=True)
        source_files = subprocess.run(
            [
                "git", "-C", REPO_ROOT, "ls-files", "-z",
                "--cached", "--others", "--exclude-standard",
            ],
            capture_output=True)
    except OSError as error:
        raise MeasurementError(f"could not identify source commit: {error}")
    commit = git.stdout.strip()
    if git.returncode != 0 or len(commit) not in (40, 64) or any(
            character not in "0123456789abcdef" for character in commit):
        raise MeasurementError("could not identify the source Git commit")
    if source_files.returncode != 0:
        raise MeasurementError("could not enumerate the source tree")

    source_digest = hashlib.sha256()
    source_count = 0
    relative_files = sorted(
        raw.decode("utf-8", "surrogateescape")
        for raw in source_files.stdout.split(b"\0") if raw)
    for relative in relative_files:
        if (not relative.endswith(_SOURCE_EXTENSIONS) and
                os.path.basename(relative) != "CMakeLists.txt"):
            continue
        path = os.path.join(REPO_ROOT, relative)
        try:
            content_sha256 = _sha256_file(path)
        except OSError as error:
            raise MeasurementError(
                f"could not identify source file {relative}: {error}")
        source_digest.update(relative.encode("utf-8", "surrogateescape"))
        source_digest.update(b"\0")
        source_digest.update(content_sha256.encode("ascii"))
        source_digest.update(b"\0")
        source_count += 1
    if source_count == 0:
        raise MeasurementError("source identity contains no source files")
    return {
        "benchmark": benchmark,
        "source": {
            "git_commit": commit,
            "state_sha256": source_digest.hexdigest(),
            "file_count": source_count,
            "generator_sha256": generator_sha256,
        },
    }


def _python_runtime_identity():
    libc_name, libc_version = platform.libc_ver()
    return {
        "implementation": sys.implementation.name,
        "version": ".".join(str(value) for value in sys.version_info[:3]),
        "cache_tag": sys.implementation.cache_tag or "",
        "libc": f"{libc_name}-{libc_version}",
    }


def capture_artifact_identity(bench, generator):
    """Capture the executable and source state before the first measurement."""
    return _artifact_identity(bench, generator)


def verify_artifact_identity(identity, bench, generator):
    """Reject publication if executable or source state changed during a run."""
    current = _artifact_identity(bench, generator)
    if current != identity:
        raise MeasurementError(
            "benchmark, generator, or source state changed during measurement")


def verify_benchmark_identity(document, bench):
    """Refuse replay with a benchmark executable different from provenance."""
    try:
        expected = document["provenance"]["benchmark"]
        actual = _benchmark_identity(bench)
    except (KeyError, MeasurementError, TypeError) as error:
        raise MeasurementError(f"could not verify benchmark identity: {error}")
    if (actual["size"] != expected["size"] or
            actual["sha256"] != expected["sha256"]):
        raise MeasurementError(
            f"benchmark identity mismatch for {actual['path']!r}; table expects "
            f"sha256={expected['sha256']}")


def derive_seed(base_seed, *domain):
    """Derive a stable, domain-separated uint64 benchmark seed."""
    if (isinstance(base_seed, bool) or not isinstance(base_seed, int) or
            not 0 <= base_seed <= 0xffffffffffffffff):
        raise ValueError("base seed must be a uint64")
    digest = hashlib.blake2b(digest_size=8, person=b"wh2-peel-v2")
    digest.update(base_seed.to_bytes(8, "little"))
    for item in domain:
        encoded = str(item).encode("utf-8")
        digest.update(len(encoded).to_bytes(4, "little"))
        digest.update(encoded)
    return int.from_bytes(digest.digest(), "little")


_STOCK_CACHE = {}


def _production_staircase_max(block_count):
    """Largest S accepted by the production mixed precode parameter domain."""
    return (
        _PRODUCTION_TOTAL_SPAN_MAX - block_count -
        _PRODUCTION_DENSE_ROWS - _PRODUCTION_MIXED_HEAVY_ROWS
    )


def _cpp_trunc_div(numerator, denominator):
    """Divide signed integers with C++'s truncation toward zero."""
    quotient = abs(numerator) // denominator
    return -quotient if numerator < 0 else quotient


def _legacy_dense_count(block_count):
    """Mirror wirehair::GetDenseCount() exactly for K in [2, 64000]."""
    if (type(block_count) is not int or
            not 2 <= block_count <= 64000):
        raise ValueError("legacy dense-count K must be in [2,64000]")

    if block_count < len(_LEGACY_TINY_DENSE_COUNTS):
        return _LEGACY_TINY_DENSE_COUNTS[block_count]
    if block_count < 2048:
        if block_count <= 500:
            low_k, high_k, low_count, high_count = 64, 500, 26, 35
        elif block_count <= 1000:
            low_k, high_k, low_count, high_count = 500, 1000, 35, 48
        else:
            low_k, high_k, low_count, high_count = 1000, 2048, 48, 62
    else:
        # Preserve native's strict ``N > mid.N`` binary search, including its
        # duplicate K=54811 point.  Equality selects the preceding segment;
        # K=54812 selects the segment beginning at the second duplicate.
        low = 0
        high = len(_LEGACY_DENSE_POINTS) - 1
        while True:
            middle = (high + low) // 2
            if middle == low:
                break
            if block_count > _LEGACY_DENSE_POINTS[middle][0]:
                low = middle
            else:
                high = middle
        low_k, low_count = _LEGACY_DENSE_POINTS[low]
        high_k, high_count = _LEGACY_DENSE_POINTS[low + 1]

    numerator = (block_count - low_k) * (high_count - low_count)
    dense_count = (
        low_count +
        _cpp_trunc_div(numerator, high_k - low_k)
    )
    # Native rounds upward to the next D such that D mod 4 == 2.
    dense_count += (2 - dense_count % 4) % 4
    return dense_count


def _dispatch_staircase_count(block_count):
    """Mirror dispatch-v1 SmallBandStaircaseCount() exactly."""
    if (type(block_count) is not int or
            block_count < 2 or block_count > 64000):
        return 0
    inherited = _legacy_dense_count(block_count)
    if block_count > 100:
        return inherited
    staircase = 2
    while 16 * (staircase + 1) ** 2 <= 25 * block_count:
        staircase += 1
    return min(staircase, inherited)


def _dispatch_source_hits(block_count):
    """Mirror dispatch-v1 CertifiedSourceHits()."""
    return 3 if block_count >= 10000 else 2


def stock_profile(bench, block_count, *, target_profile):
    """Return the exact native PMF and dispatch contract for one K."""
    if (isinstance(block_count, bool) or
            not isinstance(block_count, int) or
            not 2 <= block_count <= 64000):
        raise ValueError("native PMF K must be in [2,64000]")
    if target_profile != TARGET_PROFILE:
        raise ValueError(
            f"target_profile must be exactly {TARGET_PROFILE!r}")
    benchmark = _benchmark_identity(bench)
    key = (
        benchmark["path"], benchmark["sha256"], target_profile,
        TARGET_CONTRACT["contract_id"], TARGET_CONTRACT["contract_sha256"],
        block_count,
    )
    cached = _STOCK_CACHE.get(key)
    if cached is not None:
        return cached

    stdout = _run_checked(
        [
            bench, "peelpmf", "--N", str(block_count),
            "--target-profile", target_profile,
        ],
        isolated_codec_env())
    if _benchmark_identity(bench) != benchmark:
        raise MeasurementError(
            "benchmark changed while reading the native peel profile")
    metadata = None
    csv_header = None
    probabilities = {}
    for line in stdout.splitlines():
        if not line.strip():
            continue
        if line.startswith("# peelpmf,"):
            if metadata is not None or csv_header is not None or probabilities:
                raise MeasurementError(
                    f"misordered or duplicate peelpmf metadata "
                    f"for K={block_count}")
            fields = line.split(",")
            try:
                values = {}
                for field in fields[1:]:
                    name, value = field.split("=", 1)
                    if name in values:
                        raise ValueError("duplicate metadata field")
                    values[name] = value
                expected_fields = {
                    "N", "target_profile", "contract_id",
                    "contract_sha256", "precode_contract",
                    "packet_contract", "architecture", "degrees",
                    "staircase", "dense_rows", "heavy_rows", "source_hits",
                    "completion", "gf256_rows", "gf16_rows", "period",
                    "geometry", "residue_schedule", "residue_skew",
                    "mix_count", "target_mean", "pmf_sha256",
                    "pmf_encoding",
                }
                if set(values) != expected_fields:
                    raise ValueError("unexpected metadata fields")
                metadata = {
                    "N": int(values["N"]),
                    "target_profile": values["target_profile"],
                    "contract_id": values["contract_id"],
                    "contract_sha256": values["contract_sha256"],
                    "precode_contract": int(values["precode_contract"]),
                    "packet_contract": int(values["packet_contract"]),
                    "architecture": values["architecture"],
                    "degrees": int(values["degrees"]),
                    "staircase": int(values["staircase"]),
                    "dense_rows": int(values["dense_rows"]),
                    "heavy_rows": int(values["heavy_rows"]),
                    "source_hits": int(values["source_hits"]),
                    "completion": values["completion"],
                    "gf256_rows": int(values["gf256_rows"]),
                    "gf16_rows": int(values["gf16_rows"]),
                    "period": int(values["period"]),
                    "geometry": values["geometry"],
                    "residue_schedule": values["residue_schedule"],
                    "residue_skew": int(values["residue_skew"]),
                    "mix_count": int(values["mix_count"]),
                    "target_mean": float(values["target_mean"]),
                    "pmf_sha256": values["pmf_sha256"],
                    "pmf_encoding": values["pmf_encoding"],
                }
            except (KeyError, TypeError, ValueError):
                raise MeasurementError(
                    f"malformed peelpmf metadata for K={block_count}: {line}")
            continue
        fields = [field.strip() for field in line.split(",")]
        if fields == ["degree", "probability"]:
            if metadata is None or csv_header is not None or probabilities:
                raise MeasurementError(
                    f"misordered or duplicate peelpmf CSV header "
                    f"for K={block_count}")
            csv_header = fields
            continue
        if csv_header is None or len(fields) != len(csv_header):
            raise MeasurementError(
                f"unexpected peelpmf output for K={block_count}: {line}")
        try:
            degree = int(fields[csv_header.index("degree")])
            probability = float(fields[csv_header.index("probability")])
        except (ValueError, IndexError):
            raise MeasurementError(
                f"malformed peelpmf row for K={block_count}: {line}")
        if degree in probabilities or degree < 1 or degree > 64:
            raise MeasurementError(
                f"invalid peelpmf degree for K={block_count}: {line}")
        probabilities[degree] = probability

    if metadata is None:
        raise MeasurementError(
            f"peelpmf returned no metadata for K={block_count}")
    expected_staircase = _dispatch_staircase_count(block_count)
    expected_source_hits = _dispatch_source_hits(block_count)
    expected_target_mean = (
        block_count * min(expected_source_hits, expected_staircase) /
        expected_staircase
    )
    contract_matches = all(
        metadata[name] == value
        for name, value in TARGET_CONTRACT.items()
        if name in metadata
    )
    if (metadata["N"] != block_count or
            metadata["target_profile"] != target_profile or
            not contract_matches or metadata["degrees"] != 64 or
            metadata["staircase"] != expected_staircase or
            metadata["staircase"] >
            _production_staircase_max(block_count) or
            metadata["source_hits"] != expected_source_hits or
            not math.isfinite(metadata["target_mean"]) or
            metadata["target_mean"] <= 0.0 or
            metadata["target_mean"] != expected_target_mean or
            metadata["pmf_sha256"] != STOCK_PMF_DIGEST):
        raise MeasurementError(
            f"invalid peelpmf metadata for K={block_count}: {metadata}")
    if set(probabilities) != set(range(1, 65)):
        raise MeasurementError(
            f"peelpmf returned {len(probabilities)} of 64 degrees "
            f"for K={block_count}")
    pmf = tuple(probabilities[degree] for degree in range(1, 65))
    if (any(not math.isfinite(value) or value < 0.0 for value in pmf) or
            abs(sum(pmf) - 1.0) > 1e-12):
        raise MeasurementError(
            f"peelpmf returned a non-probability distribution "
            f"for K={block_count}")
    profile = StockProfile(
        block_count=metadata["N"],
        target_profile=metadata["target_profile"],
        contract_id=metadata["contract_id"],
        contract_sha256=metadata["contract_sha256"],
        precode_contract=metadata["precode_contract"],
        packet_contract=metadata["packet_contract"],
        architecture=metadata["architecture"],
        staircase=metadata["staircase"],
        dense_rows=metadata["dense_rows"],
        heavy_rows=metadata["heavy_rows"],
        source_hits=metadata["source_hits"],
        completion=metadata["completion"],
        gf256_rows=metadata["gf256_rows"],
        gf16_rows=metadata["gf16_rows"],
        period=metadata["period"],
        geometry=metadata["geometry"],
        residue_schedule=metadata["residue_schedule"],
        residue_skew=metadata["residue_skew"],
        mix_count=metadata["mix_count"],
        target_mean=metadata["target_mean"],
        native_pmf_sha256=metadata["pmf_sha256"],
        pmf_encoding=metadata["pmf_encoding"],
        pmf=pmf,
    )
    _STOCK_CACHE[key] = profile
    return profile


def stock_pmf(bench, block_count, *, target_profile):
    """Return degrees 1..64 of the exact PMF recovered by the native codec."""
    return list(stock_profile(
        bench, block_count, target_profile=target_profile).pmf)


def pmf_sha256(pmf):
    """Hash the exact binary64 text law passed through the benchmark hook."""
    try:
        encoded = ",".join(
            f"{float(probability):.17g}" for probability in pmf).encode("ascii")
    except (OverflowError, TypeError, ValueError):
        raise ValueError("invalid PMF")
    return hashlib.sha256(encoded).hexdigest()


def family(stock, p1, tilt, dmax, absorb):
    """Apply the native search family's coordinates to a shipped PMF."""
    try:
        coordinates = tuple(float(value)
                            for value in (p1, tilt, dmax, absorb))
    except (OverflowError, TypeError, ValueError):
        return None
    if any(not math.isfinite(value) for value in coordinates):
        return None
    p1, tilt, dmax, absorb = coordinates
    dmax = max(2, min(int(dmax), len(stock)))
    weights = [0.0] * dmax
    try:
        weights[0] = (p1 / 100.0) * stock[0]
        for degree in range(2, dmax):
            weights[degree - 1] = (
                stock[degree - 1] * (degree ** (-tilt / 100.0)))
        absorption = min(max(absorb, 0), 100) / 100.0
        weights[dmax - 1] = (
            absorption * sum(stock[dmax - 1:]) +
            (1.0 - absorption) * stock[dmax - 1]
        ) * (dmax ** (-tilt / 100.0))
    except (OverflowError, TypeError, ValueError):
        return None
    total = sum(value for value in weights if value > 0.0)
    if not math.isfinite(total) or not total > 0.0:
        return None
    return [max(value, 0.0) / total for value in weights]


def compare_probe(
        bench, block_count, trials, block_bytes, peel_weights=None,
        degree_scale=None, *, native_profile, target_profile, seed_policy,
        loss, schedule, construction_seed, loss_seed):
    """Measure one exact raw dispatch-v1 arm under split deterministic seeds."""
    for name, value in (
            ("construction_seed", construction_seed),
            ("loss_seed", loss_seed)):
        if (isinstance(value, bool) or not isinstance(value, int) or
                not 0 <= value <= 0xffffffffffffffff):
            raise ValueError(f"{name} must be a uint64")
    if target_profile != TARGET_PROFILE:
        raise ValueError(
            f"target_profile must be exactly {TARGET_PROFILE!r}")
    if seed_policy != TARGET_SEED_POLICY:
        raise ValueError(f"seed_policy must be exactly {TARGET_SEED_POLICY!r}")
    if not isinstance(schedule, str) or schedule not in TARGET_SCHEDULES:
        raise ValueError("invalid target schedule")
    if not valid_loss_rate(loss):
        raise ValueError("loss must be finite and in [0,0.99]")
    loss = float(loss)
    if (isinstance(block_count, bool) or
            not isinstance(block_count, int) or
            not 2 <= block_count <= 64000 or
            isinstance(trials, bool) or
            not isinstance(trials, int) or
            not 1 <= trials <= _COMPARE_TRIALS_MAX or
            isinstance(block_bytes, bool) or
            not isinstance(block_bytes, int) or
            not 1 <= block_bytes <= _COMPARE_BLOCK_BYTES_MAX or
            block_bytes % 2 != 0):
        raise ValueError("invalid compare K, trial count, or block size")
    if (not isinstance(native_profile, StockProfile) or
            native_profile.block_count != block_count or
            native_profile.target_profile != target_profile):
        raise ValueError("native_profile does not match compare K and target")
    _validate_native_receipt(
        native_profile.as_dict(), block_count, "compare native_profile")

    overrides = {}
    peel_spec = None
    if peel_weights is not None:
        try:
            numeric_weights = []
            for value in peel_weights:
                if isinstance(value, bool):
                    raise ValueError("Boolean peel weight")
                numeric_weights.append(float(value))
        except (OverflowError, TypeError, ValueError):
            raise ValueError("invalid peel weight vector")
        total_weight = sum(numeric_weights)
        if (not 2 <= len(numeric_weights) <= 64 or
                any(not math.isfinite(value) or value < 0.0
                    for value in numeric_weights) or
                not math.isfinite(total_weight) or not total_weight > 0.0):
            raise ValueError("invalid peel weight vector")
        peel_spec = ",".join(
            f"{value:.17g}" for value in numeric_weights)
        overrides["WIREHAIR_V2_PEEL_DEGREES"] = peel_spec
    scale_spec = "unset"
    if degree_scale is not None:
        if isinstance(degree_scale, bool):
            raise ValueError("invalid staircase degree scale")
        try:
            numeric_scale = float(degree_scale)
        except (OverflowError, TypeError, ValueError):
            raise ValueError("invalid staircase degree scale")
        if (not math.isfinite(numeric_scale) or
                not 0.0 <= numeric_scale <= 64000.0):
            raise ValueError("invalid staircase degree scale")
        scale_spec = _canonical_staircase_scale_spec(numeric_scale)
        overrides["WIREHAIR_V2_STAIRCASE_DEGREE_SCALE"] = scale_spec
    expected_pmf_digest = hashlib.sha256(
        (peel_spec if peel_spec is not None else "stock").encode(
            "ascii")).hexdigest()

    command = [
        bench, "compare",
        "--nlo", str(block_count),
        "--nhi", str(block_count),
        "--bb-list", str(block_bytes),
        "--trials", str(trials),
        "--loss", f"{loss:.17g}",
        "--loss-seed", str(loss_seed),
        "--schedule", schedule,
        # A fixed cap can silently suppress the requested row when K * bb is
        # larger than the benchmark's default 128 MiB cap.
        "--max-message-mib", "0",
        "--precode",
        "--target-profile", target_profile,
        "--seed-policy", seed_policy,
        "--construction-seed", str(construction_seed),
    ]
    stdout = _run_checked(command, isolated_codec_env(overrides))
    target_receipt = None
    compare_metadata = None
    header = None
    rows = []
    row_codecs = []
    for line in stdout.splitlines():
        if not line.strip():
            continue
        if line.startswith("# wh2_target,"):
            if (target_receipt is not None or compare_metadata is not None or
                    header is not None or row_codecs):
                raise MeasurementError(
                    f"misordered or duplicate wh2_target receipt "
                    f"for K={block_count}")
            target_receipt = {}
            for field in line[len("# wh2_target,"):].split(","):
                if "=" not in field:
                    raise MeasurementError(
                        f"malformed wh2_target receipt for K={block_count}: "
                        f"{line}")
                name, value = field.split("=", 1)
                if not name or name in target_receipt:
                    raise MeasurementError(
                        f"duplicate wh2_target field {name!r} "
                        f"for K={block_count}")
                target_receipt[name] = value
            continue
        if line.startswith("# compare:"):
            if (target_receipt is None or compare_metadata is not None or
                    header is not None or row_codecs):
                raise MeasurementError(
                    f"misordered or duplicate compare metadata "
                    f"for K={block_count}")
            compare_metadata = {}
            for field in line[len("# compare:"):].split():
                if "=" not in field:
                    raise MeasurementError(
                        f"malformed compare metadata for K={block_count}: "
                        f"{line}")
                name, value = field.split("=", 1)
                if name in compare_metadata:
                    raise MeasurementError(
                        f"duplicate compare metadata field {name!r} "
                        f"for K={block_count}")
                compare_metadata[name] = value
            continue
        if line.startswith("codec"):
            parsed_header = tuple(line.split())
            if (compare_metadata is None or header is not None or row_codecs or
                    parsed_header != COMPARE_COLUMNS):
                raise MeasurementError(
                    f"unexpected, misordered, or duplicate compare header "
                    f"for K={block_count}")
            header = list(parsed_header)
            continue
        fields = line.split()
        codec = fields[0]
        if codec not in ("baseline", "v2", "v2_target"):
            raise MeasurementError(
                f"unexpected compare output for K={block_count}: {line}")
        if (header is None or len(fields) != len(header) or
                codec in row_codecs):
            raise MeasurementError(
                f"malformed or duplicate compare row for K={block_count}: "
                f"{line}")
        row_codecs.append(codec)
        try:
            row_block_bytes = int(fields[header.index("bb")])
            row_trials = int(fields[header.index("trials")])
            row_fail = int(fields[header.index("fail")])
            row_numeric = {
                name: float(fields[header.index(name)])
                for name in (
                    "N_mean", "OH_mean", "OH_sd", "OH50", "OH95",
                    "OH99", "OH_max", "create_MBps", "encode_MBps",
                    "decode_MBps", "recover_MBps")
            }
        except (ValueError, IndexError):
            raise MeasurementError(
                f"malformed compare row for K={block_count}: {line}")
        row_n_mean = row_numeric["N_mean"]
        if (row_block_bytes != block_bytes or row_trials != trials or
                not 0 <= row_fail <= trials or
                not math.isfinite(row_n_mean) or
                row_n_mean != float(block_count) or
                any(not math.isfinite(value) or value < 0.0
                    for value in row_numeric.values()) or
                row_numeric["OH_mean"] > row_numeric["OH_max"] or
                not (row_numeric["OH50"] <= row_numeric["OH95"] <=
                     row_numeric["OH99"] <= row_numeric["OH_max"])):
            raise MeasurementError(
                f"semantically wrong compare row for K={block_count}: {line}")
        if codec == "v2_target":
            try:
                values = {
                    name: fields[header.index(name)]
                    for name in (
                        "bb", "trials", "N_mean", "fail", "OH_mean", "OH_sd",
                        "OH50", "OH95", "OH99", "OH_max", "decode_MBps")
                }
                metrics = RecoveryMetrics(
                    construction_seed=construction_seed,
                    loss_seed=loss_seed,
                    target_receipt=dict(target_receipt or {}),
                    fail=int(values["fail"]),
                    oh_mean=float(values["OH_mean"]),
                    oh_sd=float(values["OH_sd"]),
                    oh50=float(values["OH50"]),
                    oh95=float(values["OH95"]),
                    oh99=float(values["OH99"]),
                    oh_max=float(values["OH_max"]),
                    decode_mbps=float(values["decode_MBps"]),
                )
            except (KeyError, ValueError, IndexError):
                raise MeasurementError(
                    f"malformed compare row for K={block_count}: {line}")
            numeric = (
                metrics.oh_mean, metrics.oh_sd, metrics.oh50, metrics.oh95,
                metrics.oh99, metrics.oh_max, metrics.decode_mbps)
            if (not 0 <= metrics.fail <= trials or
                    any(not math.isfinite(value) or value < 0.0
                        for value in numeric) or
                    metrics.oh_mean > metrics.oh_max or
                    not (metrics.oh50 <= metrics.oh95 <= metrics.oh99 <=
                         metrics.oh_max)):
                raise MeasurementError(
                    f"invalid compare row for K={block_count}: {line}")
            rows.append(metrics)
    if target_receipt is None:
        raise MeasurementError(
            f"compare returned no wh2_target receipt for K={block_count}")
    if compare_metadata is None:
        raise MeasurementError(
            f"compare returned no metadata for K={block_count}")

    expected_target_receipt = {
        "profile": target_profile,
        "contract_id": TARGET_CONTRACT["contract_id"],
        "contract_sha256": TARGET_CONTRACT["contract_sha256"],
        "precode_contract": str(TARGET_CONTRACT["precode_contract"]),
        "packet_contract": str(TARGET_CONTRACT["packet_contract"]),
        "architecture": TARGET_CONTRACT["architecture"],
        **TARGET_MEASUREMENT_POLICY,
        "N": str(block_count),
        "bb": str(block_bytes),
        "staircase": str(native_profile.staircase),
        "dense_rows": str(native_profile.dense_rows),
        "heavy_rows": str(native_profile.heavy_rows),
        "completion": native_profile.completion,
        "gf256_rows": str(native_profile.gf256_rows),
        "gf16_rows": str(native_profile.gf16_rows),
        "period": str(native_profile.period),
        "geometry": native_profile.geometry,
        "residue_schedule": native_profile.residue_schedule,
        "residue_skew": str(native_profile.residue_skew),
        "source_hits": str(native_profile.source_hits),
        "target_mean": _target_mean_spec(
            block_count, native_profile.staircase,
            native_profile.source_hits, scale_spec),
        "mix_count": str(native_profile.mix_count),
        "packet_peel_seed": str(
            ((construction_seed & 0xffffffff) ^
             ((construction_seed >> 32) & 0xffffffff)) & 0xffffffff),
        "construction_seed": f"0x{construction_seed:x}",
        "loss_rate": f"{loss:.17g}",
        "loss_seed": f"0x{loss_seed:x}",
        "schedule": schedule,
        "pmf_sha256": expected_pmf_digest,
        "pmf_encoding": TARGET_CONTRACT["pmf_encoding"],
        "staircase_scale": scale_spec,
    }
    if target_receipt != expected_target_receipt:
        differences = sorted(
            set(target_receipt) ^ set(expected_target_receipt))
        mismatches = sorted(
            name for name in set(target_receipt) & set(expected_target_receipt)
            if target_receipt[name] != expected_target_receipt[name])
        raise MeasurementError(
            f"wh2_target receipt does not match the requested trial for "
            f"K={block_count}: fields={differences}, values={mismatches}")

    try:
        metadata_seed = int(compare_metadata["seed"], 0)
        schedule_seed = int(compare_metadata["schedule_seed"], 0)
        metadata_loss = float(compare_metadata["loss"])
        dynamic_matches = (
            compare_metadata["N"] == f"[{block_count},{block_count}]" and
            compare_metadata["trials/bb"] == str(trials) and
            compare_metadata["max_message_mib"] == "0" and
            compare_metadata["schedule"] == schedule and
            compare_metadata["loss_trace"] == "common-id-v2" and
            metadata_seed == loss_seed and schedule_seed == loss_seed and
            math.isfinite(metadata_loss) and metadata_loss == loss
        )
    except (KeyError, TypeError, ValueError):
        dynamic_matches = False
    if not dynamic_matches:
        raise MeasurementError(
            f"compare metadata does not match the requested trial for "
            f"K={block_count}")
    expected_metadata = {
        "N", "trials/bb", "loss", "seed", "max_message_mib", "schedule",
        "schedule_seed", "loss_trace",
        *PRODUCTION_COMPLETION_REGIME,
        *_COMPARE_BANNER_ARM,
    }
    if set(compare_metadata) != expected_metadata:
        raise MeasurementError(
            f"compare returned unexpected metadata fields for "
            f"K={block_count}: "
            f"{sorted(set(compare_metadata) ^ expected_metadata)}")
    completion = {
        name: compare_metadata.get(name)
        for name in PRODUCTION_COMPLETION_REGIME
    }
    if completion != PRODUCTION_COMPLETION_REGIME:
        raise MeasurementError(
            f"compare did not use production mixed completion for "
            f"K={block_count}: {completion}")
    compare_arm = {
        name: compare_metadata.get(name)
        for name in _COMPARE_BANNER_ARM
    }
    if compare_arm != _COMPARE_BANNER_ARM:
        raise MeasurementError(
            f"compare did not use the production WH2 arm for "
            f"K={block_count}: {compare_arm}")
    if row_codecs != ["baseline", "v2", "v2_target"]:
        raise MeasurementError(
            f"compare returned unexpected row sequence {row_codecs!r} "
            f"for K={block_count}")
    if len(rows) != 1:
        raise MeasurementError(
            f"compare returned {len(rows)} v2_target rows for K={block_count}")
    return rows[0]


def make_table_document(
        entries, *, generator, bench, construction_seed_base, loss_seed_base,
        target_profile, seed_policy, loss, schedule, settings,
        source_provenance=None, artifact_identity=None):
    """Wrap entries in the strict, replayable peel-table schema."""
    if (not isinstance(generator, str) or not generator or
            not isinstance(settings, dict) or
            target_profile != TARGET_PROFILE or
            seed_policy != TARGET_SEED_POLICY or
            not isinstance(schedule, str) or
            schedule not in TARGET_SCHEDULES or
            not valid_loss_rate(loss) or
            any(
                isinstance(seed, bool) or not isinstance(seed, int) or
                not 0 <= seed <= 0xffffffffffffffff
                for seed in (construction_seed_base, loss_seed_base))):
        raise ValueError(
            "invalid table generator, target, settings, loss, or seeds")
    if artifact_identity is None:
        identity = _artifact_identity(bench, generator)
    else:
        verify_artifact_identity(artifact_identity, bench, generator)
        identity = artifact_identity
    document = {
        "schema": PEEL_TABLE_SCHEMA,
        "provenance": {
            "generator": generator,
            "benchmark": identity["benchmark"],
            "source": identity["source"],
            "native_pmf": NATIVE_PMF_PROTOCOL,
            "native_compare": NATIVE_COMPARE_PROTOCOL,
            "seed_derivation": SEED_DERIVATION_PROTOCOL,
            "completion_regime": {
                "protocol": COMPLETION_REGIME_PROTOCOL,
                "settings": dict(PRODUCTION_COMPLETION_REGIME),
            },
            "target_contract": dict(TARGET_CONTRACT),
            "compare_arm": dict(PRODUCTION_COMPARE_ARM),
            "measurement_policy": {
                "target_profile": target_profile,
                "seed_policy": seed_policy,
                "loss": float(loss),
                "schedule": schedule,
            },
            "recovery_metric_scope": dict(RECOVERY_METRIC_SCOPE),
            "python_runtime": _python_runtime_identity(),
            "construction_seed_base": construction_seed_base,
            "loss_seed_base": loss_seed_base,
            "settings": settings,
        },
        "entries": {
            str(k): value for k, value in sorted(entries.items())
        },
    }
    if source_provenance is not None:
        document["source_provenance"] = source_provenance
    validate_table_document(document)
    return document


def _is_sha256(value):
    return (
        isinstance(value, str) and len(value) == 64 and
        all(character in "0123456789abcdef" for character in value)
    )


_TARGET_CONTRACT_INTEGER_FIELDS = {
    "precode_contract", "packet_contract", "dense_rows", "heavy_rows",
    "gf256_rows", "gf16_rows", "period", "residue_skew", "mix_count",
}
_NATIVE_INTEGER_FIELDS = {
    "block_count", "precode_contract", "packet_contract", "staircase",
    "dense_rows", "heavy_rows", "source_hits", "gf256_rows", "gf16_rows",
    "period", "residue_skew", "mix_count",
}
_PROXY_REGIME_INTEGER_FIELDS = {
    "solve_block_bytes", "cost_model_block_bytes", "cost_model_verified",
    "band_tracking_x", "seed_base", "period", "gf16_rows",
}


def _validate_target_contract_object(contract, label):
    """Require the canonical contract values and their canonical JSON types."""
    if (not isinstance(contract, dict) or
            set(contract) != set(TARGET_CONTRACT) or
            any(type(contract.get(name)) is not int
                for name in _TARGET_CONTRACT_INTEGER_FIELDS) or
            contract != TARGET_CONTRACT):
        raise MeasurementError(f"{label} is not the exact dispatch-v1 contract")


def _validate_proxy_measure_regime(regime, label):
    """Reject bool/float aliases in the fixed integer proxy receipt."""
    if (not isinstance(regime, dict) or
            set(regime) != set(PROXY_MEASURE_REGIME) or
            any(type(regime.get(name)) is not int
                for name in _PROXY_REGIME_INTEGER_FIELDS) or
            regime != PROXY_MEASURE_REGIME):
        raise MeasurementError(f"{label} is not the exact proxy regime")


_TARGET_RECEIPT_FIELDS = {
    "profile", "contract_id", "contract_sha256", "precode_contract",
    "packet_contract", "architecture", *TARGET_MEASUREMENT_POLICY,
    "N", "bb", "staircase", "dense_rows", "heavy_rows", "completion",
    "gf256_rows", "gf16_rows", "period", "geometry", "residue_schedule",
    "residue_skew", "source_hits", "target_mean", "mix_count",
    "packet_peel_seed",
    "construction_seed", "loss_rate", "loss_seed", "schedule",
    "pmf_sha256", "pmf_encoding", "staircase_scale",
}


def _validate_target_receipt(
        receipt, *, block_count, block_bytes, native, construction_seed,
        loss_seed, loss, schedule, pmf_digest, staircase_scale, label):
    """Validate every field in one native ``# wh2_target`` receipt."""
    if not isinstance(receipt, dict) or set(receipt) != _TARGET_RECEIPT_FIELDS:
        raise MeasurementError(f"{label} is incomplete")
    expected = {
        "profile": TARGET_PROFILE,
        "contract_id": TARGET_CONTRACT["contract_id"],
        "contract_sha256": TARGET_CONTRACT["contract_sha256"],
        "precode_contract": str(TARGET_CONTRACT["precode_contract"]),
        "packet_contract": str(TARGET_CONTRACT["packet_contract"]),
        "architecture": TARGET_CONTRACT["architecture"],
        **TARGET_MEASUREMENT_POLICY,
        "N": str(block_count),
        "bb": str(block_bytes),
        "staircase": str(native["staircase"]),
        "dense_rows": str(native["dense_rows"]),
        "heavy_rows": str(native["heavy_rows"]),
        "completion": native["completion"],
        "gf256_rows": str(native["gf256_rows"]),
        "gf16_rows": str(native["gf16_rows"]),
        "period": str(native["period"]),
        "geometry": native["geometry"],
        "residue_schedule": native["residue_schedule"],
        "residue_skew": str(native["residue_skew"]),
        "source_hits": str(native["source_hits"]),
        "target_mean": _target_mean_spec(
            block_count, native["staircase"], native["source_hits"],
            staircase_scale),
        "mix_count": str(native["mix_count"]),
        "packet_peel_seed": str(
            ((construction_seed & 0xffffffff) ^
             ((construction_seed >> 32) & 0xffffffff)) & 0xffffffff),
        "construction_seed": f"0x{construction_seed:x}",
        "loss_rate": f"{float(loss):.17g}",
        "loss_seed": f"0x{loss_seed:x}",
        "schedule": schedule,
        "pmf_sha256": pmf_digest,
        "pmf_encoding": TARGET_CONTRACT["pmf_encoding"],
        "staircase_scale": staircase_scale,
    }
    if receipt != expected:
        raise MeasurementError(f"{label} does not match its exact target arm")


def _validate_recovery_metrics(receipt, trials, label):
    if (not isinstance(receipt, dict) or
            not RECOVERY_METRIC_FIELDS.issubset(receipt)):
        raise MeasurementError(f"{label} is missing recovery metrics")
    fail = receipt["fail"]
    if (any(
            type(receipt[name]) is not int or
            not 0 <= receipt[name] <= 0xffffffffffffffff
            for name in ("construction_seed", "loss_seed")) or
            not isinstance(receipt["target_receipt"], dict) or
            set(receipt["target_receipt"]) != _TARGET_RECEIPT_FIELDS or
            type(fail) is not int or
            not 0 <= fail <= trials):
        raise MeasurementError(
            f"{label} has invalid construction/loss seeds, target, or fail")
    value_names = (
        "oh_mean", "OH_sd", "OH50", "OH95", "OH99", "OH_max",
        "decode_mbps")
    if any(type(receipt[name]) is not float for name in value_names):
        raise MeasurementError(f"{label} has non-canonical recovery metrics")
    values = [receipt[name] for name in value_names]
    if (any(not math.isfinite(value) or value < 0.0 for value in values) or
            not values[2] <= values[3] <= values[4] <= values[5] or
            values[0] > values[5]):
        raise MeasurementError(f"{label} has invalid recovery metrics")


def _canonical_goodput(metrics, block_count):
    """Recompute the exact goodput represented by recovery metrics."""
    if metrics["fail"] != 0:
        return 0.0
    return (
        metrics["decode_mbps"] * block_count /
        (block_count + metrics["oh_mean"])
    )


def _validate_search_context(
        context, block_count, search_kind, native, sampling_seed, label):
    """Validate every search-control field needed to replay one result."""
    if search_kind == "direct-real-codec":
        required = {
            "warm_start", "sampling_seed", "screen", "refine",
            "gate_trials", "gate_block_bytes", "rank_top",
        }
        warm_start = context.get("warm_start")
        valid_warm_start = warm_start is None
        if isinstance(warm_start, list) and len(warm_start) == 4 and all(
                type(value) is int for value in warm_start):
            valid_warm_start = (
                0 <= warm_start[0] <= 400 and
                -100 <= warm_start[1] <= 1600 and
                2 <= warm_start[2] <= 64 and
                0 <= warm_start[3] <= 100
            )
        integer_bounds = {
            "screen": (1, None),
            "refine": (0, None),
            "gate_trials": (1, _COMPARE_TRIALS_MAX),
            "gate_block_bytes": (2, _COMPARE_BLOCK_BYTES_MAX),
            "rank_top": (1, None),
        }
        invalid_integer = any(
            type(context.get(name)) is not int or
            context[name] < minimum or
            (maximum is not None and context[name] > maximum)
            for name, (minimum, maximum) in integer_bounds.items()
        ) if set(context) == required else True
        if (set(context) != required or invalid_integer or
                context["gate_block_bytes"] % 2 != 0 or
                context.get("sampling_seed") != sampling_seed or
                not valid_warm_start):
            raise MeasurementError(f"{label} has invalid direct-search context")
        return
    if search_kind != "unverified-proxy-funnel":
        raise MeasurementError(f"{label} has an unknown search context")
    required = {
        "proxy_cost_model", "proxy_measure_regime", "proxy_ordering",
        "search_box", "scale_centi", "warm_start", "sampling_seed",
        "screen", "refine", "finals", "screen_cells", "gate_trials",
        "gate_block_bytes", "rank_top", "threads", "batch", "cell_base",
    }
    if set(context) != required:
        raise MeasurementError(f"{label} has incomplete funnel context")
    _validate_proxy_measure_regime(
        context.get("proxy_measure_regime"),
        f"{label}.proxy_measure_regime")
    expected_scale_max = int(math.ceil(
        100.0 * min(
            float(block_count),
            max(40.0, 4.0 * float(native["target_mean"])))
    ))
    integer_bounds = {
        "screen": (1, None),
        "refine": (0, None),
        "finals": (2, None),
        "screen_cells": (1, None),
        "gate_trials": (1, _COMPARE_TRIALS_MAX),
        "gate_block_bytes": (2, _COMPARE_BLOCK_BYTES_MAX),
        "rank_top": (1, None),
        "threads": (1, 256),
        "batch": (1, None),
        "cell_base": (0, 0xffffffffffffffff),
    }
    invalid_integer = False
    for name, (minimum, maximum) in integer_bounds.items():
        value = context.get(name)
        if (type(value) is not int or
                value < minimum or
                (maximum is not None and value > maximum)):
            invalid_integer = True
            break
    warm_start = context.get("warm_start")
    valid_warm_start = warm_start is None
    if isinstance(warm_start, list) and len(warm_start) == 5 and all(
            type(value) is int for value in warm_start):
        valid_warm_start = (
            0 <= warm_start[0] <= expected_scale_max and
            0 <= warm_start[1] <= 400 and
            0 <= warm_start[2] <= 1700 and
            2 <= warm_start[3] <= 64 and
            0 <= warm_start[4] <= 100
        )
    scale_centi = context.get("scale_centi")
    valid_scale_centi = (
        isinstance(scale_centi, list) and
        len(scale_centi) == 2 and
        all(type(value) is int for value in scale_centi)
    )
    if (invalid_integer or
            context["cell_base"] + context["screen_cells"] > 0x100000000 or
            context["gate_block_bytes"] % 2 != 0 or
            context.get("proxy_cost_model") != PROXY_COST_MODEL or
            context.get("proxy_ordering") != PROXY_ORDERING_PROTOCOL or
            context.get("search_box") != SEARCH_BOX_PROTOCOL or
            not valid_scale_centi or
            scale_centi != [0, expected_scale_max] or
            context.get("sampling_seed") != sampling_seed or
            not valid_warm_start):
        raise MeasurementError(f"{label} has invalid funnel context")


def _validate_search_receipt(
        receipt, block_count, expected_construction_seed_base,
        expected_loss_seed_base, expected_seed_domain, expected_search_kind,
        native, measurement_policy, label):
    required = {
        "mode", "goodput", "trials", "block_bytes", "search_kind",
        "construction_seed_base", "loss_seed_base", "seed_domain",
        "coordinates", "peel_pmf", "peel_pmf_sha256", "shipped_control",
        "shipped_goodput", "context",
    }
    if (not isinstance(receipt, dict) or
            set(receipt) != required | RECOVERY_METRIC_FIELDS):
        raise MeasurementError(f"{label} is incomplete")
    trials = receipt["trials"]
    block_bytes = receipt["block_bytes"]
    if (receipt["mode"] not in ("trained", "scale-only", "shipped") or
            type(trials) is not int or
            not 1 <= trials <= _COMPARE_TRIALS_MAX or
            type(block_bytes) is not int or
            not 1 <= block_bytes <= _COMPARE_BLOCK_BYTES_MAX or
            block_bytes % 2 != 0 or
            receipt["search_kind"] != expected_search_kind or
            any(
                type(receipt[name]) is not int or
                receipt[name] != expected or
                not 0 <= receipt[name] <= 0xffffffffffffffff
                for name, expected in (
                    ("construction_seed_base",
                     expected_construction_seed_base),
                    ("loss_seed_base", expected_loss_seed_base))) or
            receipt["seed_domain"] != expected_seed_domain or
            not isinstance(receipt["coordinates"], dict) or
            not isinstance(receipt["context"], dict)):
        raise MeasurementError(f"{label} has invalid search metadata")
    expected_sampling_seed = derive_seed(
        expected_construction_seed_base, expected_seed_domain, block_count,
        "sampling")
    sampling_seed = receipt["context"].get("sampling_seed")
    if (type(sampling_seed) is not int or
            sampling_seed != expected_sampling_seed):
        raise MeasurementError(
            f"{label} has an invalid derived sampling seed")
    _validate_search_context(
        receipt["context"], block_count, receipt["search_kind"], native,
        expected_sampling_seed, f"{label}.context")
    coordinate_names = ("scale", "p1", "tilt", "dmax", "absorb")
    coordinates = receipt["coordinates"]
    if set(coordinates) != set(coordinate_names):
        raise MeasurementError(f"{label} has invalid search coordinates")
    scale_value = coordinates["scale"]
    integer_values = [
        coordinates[name] for name in ("p1", "tilt", "dmax", "absorb")
    ]
    if (type(scale_value) is not float or
            any(type(value) is not int for value in integer_values)):
        raise MeasurementError(f"{label} has non-numeric search coordinates")
    scale = scale_value
    p1 = coordinates["p1"]
    tilt = coordinates["tilt"]
    dmax = coordinates["dmax"]
    absorb = coordinates["absorb"]
    if (not math.isfinite(scale) or
            (scale != -1.0 and not 0.0 <= scale <= 64000.0) or
            not 0 <= p1 <= 400 or
            not -100 <= tilt <= 1600 or
            not 2 <= dmax <= 64 or
            not 0 <= absorb <= 100 or
            (expected_search_kind == "direct-real-codec" and
             scale != -1.0) or
            (expected_search_kind == "unverified-proxy-funnel" and
             receipt["mode"] in ("trained", "scale-only") and
             (not receipt["context"]["scale_centi"][0] <=
              int(round(scale * 100.0)) <=
              receipt["context"]["scale_centi"][1] or
              scale != int(round(scale * 100.0)) / 100.0))):
        raise MeasurementError(f"{label} has invalid search coordinates")
    peel_pmf_digest = _validate_pmf(
        receipt["peel_pmf"], f"{label}.peel_pmf")
    if receipt["peel_pmf_sha256"] != peel_pmf_digest:
        raise MeasurementError(f"{label} has an invalid peel PMF digest")
    if receipt["mode"] in ("shipped", "scale-only"):
        expected_pmf = list(native["pmf"])
        expected_scale = -1.0 if receipt["mode"] == "shipped" else scale
        if (scale != expected_scale or p1 != 100 or tilt != 0 or
                dmax != 64 or absorb != 100):
            raise MeasurementError(
                f"{label} has invalid stock-control coordinates")
        if receipt["mode"] == "scale-only" and scale < 0.0:
            raise MeasurementError(
                f"{label} has an unset scale-only coordinate")
    else:
        expected_pmf = family(native["pmf"], p1, tilt, dmax, absorb)
    if expected_pmf is None or receipt["peel_pmf"] != expected_pmf:
        raise MeasurementError(
            f"{label} peel PMF does not match its coordinates")
    # Funnel canonicalizes a native-shaped PMF with an active staircase scale
    # to the stock-hook scale-only arm. A trained receipt for the same PMF is
    # therefore not an artifact the current generators can produce.
    if (receipt["mode"] == "trained" and
            receipt["peel_pmf"] == list(native["pmf"])):
        raise MeasurementError(
            f"{label} labels the stock PMF as a trained arm")
    if type(receipt["goodput"]) is not float:
        raise MeasurementError(f"{label} has invalid goodput")
    goodput = receipt["goodput"]
    if not math.isfinite(goodput) or goodput < 0.0:
        raise MeasurementError(f"{label} has invalid goodput")
    _validate_recovery_metrics(receipt, trials, label)
    if receipt["fail"] != 0:
        raise MeasurementError(f"{label} selected a non-decoding search arm")
    _validate_recovery_metrics(
        receipt["shipped_control"], trials, f"{label}.shipped_control")
    if set(receipt["shipped_control"]) != RECOVERY_METRIC_FIELDS:
        raise MeasurementError(
            f"{label} shipped control has unexpected fields")
    expected_construction_seed = derive_seed(
        expected_construction_seed_base, receipt["seed_domain"], block_count,
        "rank", trials, block_bytes, "construction")
    expected_loss_seed = derive_seed(
        expected_loss_seed_base, receipt["seed_domain"], block_count, "rank",
        trials, block_bytes, "loss")
    if (receipt["construction_seed"] != expected_construction_seed or
            receipt["loss_seed"] != expected_loss_seed):
        raise MeasurementError(
            f"{label} has invalid derived construction/loss seeds")
    if (receipt["shipped_control"]["construction_seed"] !=
            expected_construction_seed or
            receipt["shipped_control"]["loss_seed"] != expected_loss_seed):
        raise MeasurementError(
            f"{label} control did not use the paired construction/loss seeds")
    scale_spec = (
        "unset" if scale == -1.0 else
        _canonical_staircase_scale_spec(scale))
    _validate_target_receipt(
        receipt["target_receipt"],
        block_count=block_count,
        block_bytes=block_bytes,
        native=native,
        construction_seed=expected_construction_seed,
        loss_seed=expected_loss_seed,
        loss=measurement_policy["loss"],
        schedule=measurement_policy["schedule"],
        pmf_digest=(
            STOCK_PMF_DIGEST
            if receipt["mode"] in ("shipped", "scale-only") else
            pmf_sha256(receipt["peel_pmf"])),
        staircase_scale=scale_spec,
        label=f"{label}.target_receipt",
    )
    _validate_target_receipt(
        receipt["shipped_control"]["target_receipt"],
        block_count=block_count,
        block_bytes=block_bytes,
        native=native,
        construction_seed=expected_construction_seed,
        loss_seed=expected_loss_seed,
        loss=measurement_policy["loss"],
        schedule=measurement_policy["schedule"],
        pmf_digest=STOCK_PMF_DIGEST,
        staircase_scale="unset",
        label=f"{label}.shipped_control.target_receipt",
    )
    candidate_target = dict(receipt["target_receipt"])
    control_target = dict(receipt["shipped_control"]["target_receipt"])
    for field in ("pmf_sha256", "staircase_scale", "target_mean"):
        candidate_target.pop(field)
        control_target.pop(field)
    if candidate_target != control_target:
        raise MeasurementError(
            f"{label} candidate/control changed more than the PMF or scale")
    expected_goodput = _canonical_goodput(receipt, block_count)
    if goodput != expected_goodput:
        raise MeasurementError(f"{label} has an inconsistent goodput")
    shipped = receipt["shipped_control"]
    expected_shipped_goodput = _canonical_goodput(shipped, block_count)
    if (type(receipt["shipped_goodput"]) is not float or
            receipt["shipped_goodput"] != expected_shipped_goodput):
        raise MeasurementError(
            f"{label} has an inconsistent shipped-control goodput")
    if (receipt["mode"] in ("trained", "scale-only") and
            not expected_goodput > expected_shipped_goodput):
        raise MeasurementError(
            f"{label} selected a candidate that did not beat shipped")
    if receipt["mode"] == "shipped" and any(
            receipt[name] != shipped[name]
            for name in (
                "construction_seed", "loss_seed", "target_receipt", "fail",
                "oh_mean", "OH_sd", "OH50", "OH95", "OH99", "OH_max",
                "decode_mbps")):
        raise MeasurementError(
            f"{label} shipped selection contradicts its control metrics")


def _validate_native_receipt(receipt, block_count, label):
    required = {
        "block_count", "target_profile", "contract_id", "contract_sha256",
        "precode_contract", "packet_contract", "architecture", "staircase",
        "dense_rows", "heavy_rows", "source_hits", "completion",
        "gf256_rows", "gf16_rows", "period", "geometry",
        "residue_schedule", "residue_skew", "mix_count", "target_mean",
        "native_pmf_sha256", "pmf_encoding", "pmf_sha256", "pmf",
    }
    if not isinstance(receipt, dict) or set(receipt) != required:
        raise MeasurementError(f"{label} is incomplete")
    if type(block_count) is not int or not 2 <= block_count <= 64000:
        raise MeasurementError(f"{label} has invalid block count")
    if any(type(receipt.get(name)) is not int
           for name in _NATIVE_INTEGER_FIELDS):
        raise MeasurementError(f"{label} has non-canonical integer metadata")
    if type(receipt["target_mean"]) is not float:
        raise MeasurementError(f"{label} has invalid target mean")
    if (not isinstance(receipt["pmf"], list) or
            any(type(value) is not float for value in receipt["pmf"])):
        raise MeasurementError(f"{label} has a non-canonical native PMF")
    target_mean = float(receipt["target_mean"])
    pmf_digest = _validate_pmf(receipt["pmf"], f"{label}.pmf")
    expected_staircase = _dispatch_staircase_count(block_count)
    expected_source_hits = _dispatch_source_hits(block_count)
    expected_target_mean = (
        block_count * min(expected_source_hits, expected_staircase) /
        expected_staircase
    )
    contract_matches = all(
        receipt.get(name) == value
        for name, value in TARGET_CONTRACT.items())
    if (receipt["block_count"] != block_count or
            type(receipt["block_count"]) is not int or
            not contract_matches or
            type(receipt["staircase"]) is not int or
            receipt["staircase"] != expected_staircase or
            receipt["staircase"] >
            _production_staircase_max(block_count) or
            type(receipt["source_hits"]) is not int or
            receipt["source_hits"] != expected_source_hits or
            not math.isfinite(target_mean) or target_mean <= 0.0 or
            target_mean != expected_target_mean or
            len(receipt["pmf"]) != 64 or
            receipt["native_pmf_sha256"] != STOCK_PMF_DIGEST or
            receipt["pmf_sha256"] != pmf_digest):
        raise MeasurementError(f"{label} has invalid native metadata")


def _validate_pmf(pmf, label):
    if not isinstance(pmf, list) or not 2 <= len(pmf) <= 64:
        raise MeasurementError(f"{label} must contain 2..64 probabilities")
    if any(type(value) is not float for value in pmf):
        raise MeasurementError(
            f"{label} contains non-canonical probabilities")
    values = list(pmf)
    if (any(not math.isfinite(value) or value < 0.0 for value in values) or
            abs(sum(values) - 1.0) > 1e-12):
        raise MeasurementError(f"{label} is not a probability distribution")
    return pmf_sha256(values)


def _validate_validation_receipt(
        receipt, block_count, construction_seed_base, loss_seed_base, native,
        measurement_policy, label):
    required = {
        "verdict", "margin_percent", "trials", "block_bytes",
        "scale", "trained_pmf_sha256", "trained_goodput",
        "shipped_goodput", "trained", "shipped",
    }
    if not isinstance(receipt, dict) or set(receipt) != required:
        raise MeasurementError(f"{label} is incomplete")
    trials = receipt["trials"]
    block_bytes = receipt["block_bytes"]
    if (receipt["verdict"] not in ("keep", "control") or
            type(trials) is not int or
            not 1 <= trials <= _COMPARE_TRIALS_MAX or
            type(block_bytes) is not int or
            not 1 <= block_bytes <= _COMPARE_BLOCK_BYTES_MAX or
            block_bytes % 2 != 0):
        raise MeasurementError(f"{label} has invalid validation metadata")
    value_names = (
        "margin_percent", "scale", "trained_goodput", "shipped_goodput")
    if any(type(receipt[name]) is not float for name in value_names):
        raise MeasurementError(
            f"{label} has non-canonical validation values")
    values = [receipt[name] for name in value_names]
    if (any(not math.isfinite(value) for value in values) or
            values[0] < 0.0 or values[1] < -1.0 or
            any(value < 0.0 for value in values[2:]) or
            not _is_sha256(receipt.get("trained_pmf_sha256"))):
        raise MeasurementError(f"{label} has invalid validation values")
    _validate_recovery_metrics(receipt["trained"], trials, f"{label}.trained")
    _validate_recovery_metrics(receipt["shipped"], trials, f"{label}.shipped")
    if (set(receipt["trained"]) != RECOVERY_METRIC_FIELDS or
            set(receipt["shipped"]) != RECOVERY_METRIC_FIELDS):
        raise MeasurementError(f"{label} has unexpected arm fields")
    if (receipt["trained"]["construction_seed"] !=
            receipt["shipped"]["construction_seed"] or
            receipt["trained"]["loss_seed"] !=
            receipt["shipped"]["loss_seed"]):
        raise MeasurementError(
            f"{label} did not use paired construction/loss seeds")
    if receipt["trained"]["fail"] > 0 and receipt["shipped"]["fail"] > 0:
        raise MeasurementError(f"{label} records two failed recovery arms")
    expected_construction_seed = derive_seed(
        construction_seed_base, "validation", block_count, trials,
        block_bytes, "construction")
    expected_loss_seed = derive_seed(
        loss_seed_base, "validation", block_count, trials, block_bytes, "loss")
    if (receipt["trained"]["construction_seed"] !=
            expected_construction_seed or
            receipt["trained"]["loss_seed"] != expected_loss_seed):
        raise MeasurementError(
            f"{label} has invalid derived construction/loss seeds")
    scale = receipt.get("scale")
    if (type(scale) is not float or
            not math.isfinite(scale) or
            (scale != -1.0 and not 0.0 <= scale <= 64000.0)):
        raise MeasurementError(f"{label} has invalid paired scale")
    for arm in ("trained", "shipped"):
        pmf_digest = (
            receipt[f"{arm}_pmf_sha256"]
            if arm == "trained" else STOCK_PMF_DIGEST)
        staircase_scale = (
            "unset" if arm == "shipped" or scale == -1.0
            else _canonical_staircase_scale_spec(scale))
        _validate_target_receipt(
            receipt[arm]["target_receipt"],
            block_count=block_count,
            block_bytes=block_bytes,
            native=native,
            construction_seed=expected_construction_seed,
            loss_seed=expected_loss_seed,
            loss=measurement_policy["loss"],
            schedule=measurement_policy["schedule"],
            pmf_digest=pmf_digest,
            staircase_scale=staircase_scale,
            label=f"{label}.{arm}.target_receipt",
        )
    trained_target = dict(receipt["trained"]["target_receipt"])
    shipped_target = dict(receipt["shipped"]["target_receipt"])
    for field in ("pmf_sha256", "staircase_scale", "target_mean"):
        trained_target.pop(field)
        shipped_target.pop(field)
    if trained_target != shipped_target:
        raise MeasurementError(
            f"{label} candidate/control changed more than the PMF or scale")
    canonical_goodputs = {}
    for arm in ("trained", "shipped"):
        metrics = receipt[arm]
        expected_goodput = _canonical_goodput(metrics, block_count)
        canonical_goodputs[arm] = expected_goodput
        if receipt[f"{arm}_goodput"] != expected_goodput:
            raise MeasurementError(
                f"{label} has inconsistent {arm} goodput")
    if (receipt["verdict"] == "keep" and
            not canonical_goodputs["trained"] >
            canonical_goodputs["shipped"] *
            (1.0 + receipt["margin_percent"] / 100.0)):
        raise MeasurementError(
            f"{label} keep verdict does not clear its margin")


_DIRECT_SETTINGS_FIELDS = {
    "kmin", "kmax", "screen", "refine", "gate_trials",
    "gate_block_bytes", "rank_trials", "rank_block_bytes", "rank_top",
    "target_profile", "seed_policy", "loss", "schedule",
}
_SWEEP_SETTINGS_FIELDS = {
    "proxy_k_ladder", "real_trials_override", "target_profile",
    "seed_policy", "loss", "schedule", "proxy_cost_model",
    "proxy_measure_regime", "proxy_ordering",
    "allow_unverified_cost_model", "search_box",
}
_ENTRY_BASE_FIELDS = {
    "K", "scale", "p1", "tilt", "dmax", "absorb",
    *RECOVERY_METRIC_FIELDS,
    "goodput", "native", "peel_pmf", "search_receipt",
}
_DIRECT_ENTRY_DIAGNOSTIC_FIELDS = {"seconds", "probes"}
_SWEEP_ENTRY_DIAGNOSTIC_FIELDS = {
    "S", "source_hits", "target_mean", "seconds", "screen",
    "screen_cells", "finals", "rejected", "real_trials",
}
_VALIDATION_ENTRY_FIELDS = {
    "validation_receipt", "verified_mbps", "verified_oh",
    "shipped_mbps", "gain_pct",
}


def _settings_match_measurement_policy(settings, measurement_policy):
    return (
        settings.get("target_profile") ==
            measurement_policy["target_profile"] and
        settings.get("seed_policy") == measurement_policy["seed_policy"] and
        type(settings.get("loss")) is float and
        settings["loss"] == measurement_policy["loss"] and
        settings.get("schedule") == measurement_policy["schedule"]
    )


def _validate_search_generator_settings(
        generator, settings, measurement_policy, label):
    """Validate the exact settings schema emitted by each table generator."""
    if not isinstance(settings, dict):
        raise MeasurementError(f"{label} is not an object")
    if generator == "tools/peel_direct.py":
        integer_fields = {
            "kmin", "kmax", "screen", "refine", "gate_trials",
            "gate_block_bytes", "rank_trials", "rank_block_bytes", "rank_top",
        }
        if (set(settings) != _DIRECT_SETTINGS_FIELDS or
                any(type(settings.get(name)) is not int
                    for name in integer_fields) or
                not 2 <= settings["kmin"] <= settings["kmax"] <= 64000 or
                settings["screen"] < 1 or settings["refine"] < 0 or
                not 1 <= settings["gate_trials"] <= _COMPARE_TRIALS_MAX or
                not 1 <= settings["rank_trials"] <= _COMPARE_TRIALS_MAX or
                not 2 <= settings["gate_block_bytes"] <=
                    _COMPARE_BLOCK_BYTES_MAX or
                settings["gate_block_bytes"] % 2 != 0 or
                not 2 <= settings["rank_block_bytes"] <=
                    _COMPARE_BLOCK_BYTES_MAX or
                settings["rank_block_bytes"] % 2 != 0 or
                settings["rank_top"] < 1 or
                not _settings_match_measurement_policy(
                    settings, measurement_policy)):
            raise MeasurementError(
                f"{label} is not an exact direct-search settings receipt")
        return
    if generator != "tools/peel_sweep.py":
        raise MeasurementError(
            f"peel table has unsupported search generator {generator!r}")
    ladder = settings.get("proxy_k_ladder")
    if (set(settings) != _SWEEP_SETTINGS_FIELDS or
            not isinstance(ladder, list) or not ladder or
            any(type(block_count) is not int or
                block_count not in PROXY_K_LADDER
                for block_count in ladder) or
            ladder != sorted(set(ladder)) or
            ladder != [
                block_count for block_count in PROXY_K_LADDER
                if block_count <= ladder[-1]
            ] or
            type(settings.get("real_trials_override")) is not int or
            not 0 <= settings["real_trials_override"] <= _COMPARE_TRIALS_MAX or
            settings.get("proxy_cost_model") != PROXY_COST_MODEL or
            settings.get("proxy_ordering") != PROXY_ORDERING_PROTOCOL or
            settings.get("allow_unverified_cost_model") is not True or
            settings.get("search_box") != SEARCH_BOX_PROTOCOL or
            not _settings_match_measurement_policy(
                settings, measurement_policy)):
        raise MeasurementError(
            f"{label} is not an exact proxy-sweep settings receipt")
    _validate_proxy_measure_regime(
        settings.get("proxy_measure_regime"),
        f"{label}.proxy_measure_regime")


def _sweep_budget(block_count):
    if block_count <= 1024:
        return 3000, 400, 16, 1000
    if block_count <= 8192:
        return 1200, 250, 12, 1000
    if block_count <= 32768:
        return 500, 150, 8, 600
    return 250, 100, 6, 400


def _sweep_real_trials(block_count, override):
    if override:
        return override
    if block_count <= 512:
        return 800
    if block_count <= 4096:
        return 400
    if block_count <= 16384:
        return 200
    return 100


def _validate_entry_generator_settings(
        generator, settings, receipt, block_count, label):
    """Cross-check one selected receipt against its generator invocation."""
    context = receipt["context"]
    if generator == "tools/peel_direct.py":
        expected_context = {
            "screen": settings["screen"],
            "refine": settings["refine"],
            "gate_trials": settings["gate_trials"],
            "gate_block_bytes": settings["gate_block_bytes"],
            "rank_top": settings["rank_top"],
        }
        if (receipt["trials"] != settings["rank_trials"] or
                receipt["block_bytes"] != settings["rank_block_bytes"] or
                any(context.get(name) != value
                    for name, value in expected_context.items())):
            raise MeasurementError(
                f"{label} contradicts its direct-search settings")
        return
    screen, refine, finals, screen_cells = _sweep_budget(block_count)
    expected_context = {
        "proxy_cost_model": settings["proxy_cost_model"],
        "proxy_measure_regime": settings["proxy_measure_regime"],
        "proxy_ordering": settings["proxy_ordering"],
        "search_box": settings["search_box"],
        "screen": screen,
        "refine": refine,
        "finals": finals,
        "screen_cells": screen_cells,
        "gate_trials": 25,
        "gate_block_bytes": 64,
        "rank_top": 3,
        "threads": 64,
        "batch": 60,
        "cell_base": 900_000_000,
    }
    if (receipt["trials"] != _sweep_real_trials(
            block_count, settings["real_trials_override"]) or
            receipt["block_bytes"] != 4096 or
            any(context.get(name) != value
                for name, value in expected_context.items())):
        raise MeasurementError(
            f"{label} contradicts its proxy-sweep settings")


def _validate_entry_field_schema(
        entry, source_generator, native, search_receipt, *,
        validated, selected_shipped, label):
    """Require the exact top-level fields emitted for this entry state."""
    diagnostics = (
        _DIRECT_ENTRY_DIAGNOSTIC_FIELDS
        if source_generator == "tools/peel_direct.py" else
        _SWEEP_ENTRY_DIAGNOSTIC_FIELDS
    )
    expected = set(_ENTRY_BASE_FIELDS)
    if validated:
        expected.update(_VALIDATION_ENTRY_FIELDS)
        if selected_shipped:
            expected.update(
                {"reverted_to_shipped", "search_would_have_lost_pct"})
        else:
            expected.update(diagnostics)
    else:
        expected.update(diagnostics)
        if search_receipt["mode"] == "shipped":
            expected.add("reverted_to_shipped")
    if set(entry) != expected:
        differences = sorted(set(entry) ^ expected)
        raise MeasurementError(
            f"{label} has unexpected/missing top-level fields {differences}")
    if ("reverted_to_shipped" in expected and
            entry["reverted_to_shipped"] is not True):
        raise MeasurementError(
            f"{label} has an invalid shipped-selection marker")

    # A validated control deliberately drops generator diagnostics when it
    # reconstructs the selected shipped arm.  All other states retain the
    # diagnostics emitted by the source generator.
    if validated and selected_shipped:
        return
    seconds = entry["seconds"]
    if (type(seconds) is not float or
            not math.isfinite(seconds) or seconds < 0.0):
        raise MeasurementError(f"{label} has invalid elapsed-time diagnostics")
    if source_generator == "tools/peel_direct.py":
        if type(entry["probes"]) is not int or entry["probes"] < 0:
            raise MeasurementError(f"{label} has invalid probe diagnostics")
        return
    context = search_receipt["context"]
    if (type(entry["S"]) is not int or
            entry["S"] != native["staircase"] or
            type(entry["source_hits"]) is not int or
            entry["source_hits"] != native["source_hits"] or
            type(entry["target_mean"]) is not float or
            entry["target_mean"] != native["target_mean"] or
            any(type(entry[name]) is not int
                for name in (
                    "screen", "screen_cells", "finals", "rejected",
                    "real_trials")) or
            entry["screen"] != context["screen"] or
            entry["screen_cells"] != context["screen_cells"] or
            entry["finals"] != context["finals"] or
            not 0 <= entry["rejected"] <= 2 * entry["finals"] or
            entry["real_trials"] != search_receipt["trials"]):
        raise MeasurementError(f"{label} has invalid proxy-sweep diagnostics")


def _expected_source_k(
        generator, settings, source_selection, validation_settings):
    """Return the complete and selected K sets implied by provenance."""
    if generator == "tools/peel_direct.py":
        complete = list(range(settings["kmin"], settings["kmax"] + 1))
    else:
        complete = list(settings["proxy_k_ladder"])
    if source_selection is None:
        return complete, complete
    selected = [
        block_count for block_count in complete
        if block_count <= validation_settings["kmax"]
    ]
    return complete, selected


def _validate_table_document_impl(document):
    """Implement strict table validation behind the public numeric boundary."""
    if not isinstance(document, dict):
        raise MeasurementError("peel table root must be an object")
    if document.get("schema") != PEEL_TABLE_SCHEMA:
        raise MeasurementError(
            f"unsupported or legacy peel table schema "
            f"{document.get('schema')!r}; expected {PEEL_TABLE_SCHEMA!r}")
    generator = (
        document.get("provenance", {}).get("generator")
        if isinstance(document.get("provenance"), dict) else None)
    expected_document_fields = {"schema", "provenance", "entries"}
    if generator == "tools/peel_validate.py":
        expected_document_fields.add("source_provenance")
    if set(document) != expected_document_fields:
        raise MeasurementError("peel table has unexpected root fields")
    provenance = document.get("provenance")
    if not isinstance(provenance, dict):
        raise MeasurementError("peel table has no provenance object")
    required = {
        "generator", "benchmark", "native_pmf", "native_compare",
        "seed_derivation", "completion_regime", "target_contract",
        "compare_arm", "measurement_policy", "construction_seed_base",
        "loss_seed_base", "settings", "source", "recovery_metric_scope",
        "python_runtime",
    }
    if set(provenance) != required:
        differences = sorted(required ^ set(provenance))
        raise MeasurementError(
            f"peel table provenance has unexpected/missing {differences}")
    if provenance["native_pmf"] != NATIVE_PMF_PROTOCOL:
        raise MeasurementError("peel table used an unsupported PMF source")
    if provenance["native_compare"] != NATIVE_COMPARE_PROTOCOL:
        raise MeasurementError("peel table used an unsupported compare protocol")
    if provenance["seed_derivation"] != SEED_DERIVATION_PROTOCOL:
        raise MeasurementError(
            "peel table used an unsupported seed-derivation protocol")
    completion_regime = provenance["completion_regime"]
    if (not isinstance(completion_regime, dict) or
            set(completion_regime) != {"protocol", "settings"} or
            completion_regime.get("protocol") != COMPLETION_REGIME_PROTOCOL or
            completion_regime.get("settings") !=
            PRODUCTION_COMPLETION_REGIME):
        raise MeasurementError(
            "peel table did not use the production mixed-completion regime")
    if provenance["compare_arm"] != PRODUCTION_COMPARE_ARM:
        raise MeasurementError(
            "peel table did not use the production WH2 compare arm")
    _validate_target_contract_object(
        provenance["target_contract"], "peel table target contract")
    measurement_policy = provenance["measurement_policy"]
    if (not isinstance(measurement_policy, dict) or
            set(measurement_policy) !=
            {"target_profile", "seed_policy", "loss", "schedule"} or
            measurement_policy.get("target_profile") != TARGET_PROFILE or
            measurement_policy.get("seed_policy") != TARGET_SEED_POLICY or
            not isinstance(measurement_policy.get("schedule"), str) or
            measurement_policy.get("schedule") not in TARGET_SCHEDULES or
            type(measurement_policy.get("loss")) is not float or
            not valid_loss_rate(measurement_policy.get("loss"))):
        raise MeasurementError(
            "peel table has an invalid target measurement policy")
    if provenance["recovery_metric_scope"] != RECOVERY_METRIC_SCOPE:
        raise MeasurementError(
            "peel table has an unsupported recovery-metric scope")
    python_runtime = provenance["python_runtime"]
    if (not isinstance(python_runtime, dict) or
            set(python_runtime) !=
            {"implementation", "version", "cache_tag", "libc"} or
            any(not isinstance(value, str)
                for value in python_runtime.values()) or
            not python_runtime["implementation"] or
            not python_runtime["version"]):
        raise MeasurementError("peel table has invalid Python runtime identity")
    benchmark = provenance["benchmark"]
    source = provenance["source"]
    if (not isinstance(provenance["generator"], str) or
            not provenance["generator"] or
            not isinstance(benchmark, dict) or
            set(benchmark) != {"path", "sha256", "size"} or
            not isinstance(benchmark.get("path"), str) or
            not benchmark["path"] or not _is_sha256(benchmark.get("sha256")) or
            type(benchmark.get("size")) is not int or
            benchmark["size"] < 1 or
            not isinstance(source, dict) or
            set(source) != {
                "git_commit", "state_sha256", "file_count",
                "generator_sha256"} or
            not isinstance(source.get("git_commit"), str) or
            len(source["git_commit"]) not in (40, 64) or
            any(character not in "0123456789abcdef"
                for character in source["git_commit"]) or
            not _is_sha256(source.get("state_sha256")) or
            not _is_sha256(source.get("generator_sha256")) or
            type(source.get("file_count")) is not int or
            source["file_count"] < 1 or
            not isinstance(provenance["settings"], dict) or
            any(
                type(provenance[name]) is not int or
                not 0 <= provenance[name] <= 0xffffffffffffffff
                for name in ("construction_seed_base", "loss_seed_base"))):
        raise MeasurementError("peel table has invalid provenance fields")
    expected_search_construction_seed_base = (
        provenance["construction_seed_base"])
    expected_search_loss_seed_base = provenance["loss_seed_base"]
    source_generator = provenance["generator"]
    source_settings = provenance["settings"]
    validation_settings = None
    source_selection = None
    if provenance["generator"] == "tools/peel_validate.py":
        source_provenance = document.get("source_provenance")
        if (not isinstance(source_provenance, dict) or
                set(source_provenance) != {
                    "schema", "document_sha256", "provenance",
                    "entry_count", "selected_entry_count", "selected_K"} or
                source_provenance.get("schema") != PEEL_TABLE_SCHEMA or
                not _is_sha256(source_provenance.get("document_sha256")) or
                not isinstance(source_provenance.get("provenance"), dict) or
                type(source_provenance.get("entry_count")) is not int or
                source_provenance["entry_count"] < 1 or
                type(source_provenance.get("selected_entry_count")) is not
                    int or
                source_provenance["selected_entry_count"] < 1 or
                not isinstance(source_provenance.get("selected_K"), list) or
                any(type(k) is not int or
                    not 2 <= k <= 64000
                    for k in source_provenance["selected_K"]) or
                source_provenance["selected_K"] !=
                sorted(set(source_provenance["selected_K"]))):
            raise MeasurementError(
                "validated peel table has an invalid source-table receipt")
        source_selection = source_provenance
        source_receipt = source_provenance["provenance"]
        source_benchmark = source_receipt.get("benchmark")
        source_identity = source_receipt.get("source")
        source_runtime = source_receipt.get("python_runtime")
        _validate_target_contract_object(
            source_receipt.get("target_contract"),
            "validated peel table source target contract")
        if (set(source_receipt) != required or
                not isinstance(source_receipt.get("generator"), str) or
                not source_receipt["generator"] or
                not isinstance(source_benchmark, dict) or
                set(source_benchmark) != {"path", "sha256", "size"} or
                not isinstance(source_benchmark.get("path"), str) or
                not source_benchmark["path"] or
                not _is_sha256(source_benchmark.get("sha256")) or
                type(source_benchmark.get("size")) is not int or
                source_benchmark["size"] < 1 or
                not isinstance(source_identity, dict) or
                set(source_identity) != {
                    "git_commit", "state_sha256", "file_count",
                    "generator_sha256"} or
                not isinstance(source_identity.get("git_commit"), str) or
                len(source_identity["git_commit"]) not in (40, 64) or
                any(character not in "0123456789abcdef"
                    for character in source_identity["git_commit"]) or
                not _is_sha256(source_identity.get("state_sha256")) or
                not _is_sha256(source_identity.get("generator_sha256")) or
                type(source_identity.get("file_count")) is not int or
                source_identity["file_count"] < 1 or
                not isinstance(source_receipt.get("settings"), dict) or
                not isinstance(source_runtime, dict) or
                set(source_runtime) !=
                {"implementation", "version", "cache_tag", "libc"} or
                any(not isinstance(value, str)
                    for value in source_runtime.values()) or
                not source_runtime["implementation"] or
                not source_runtime["version"] or
                source_receipt.get("native_pmf") != NATIVE_PMF_PROTOCOL or
                source_receipt.get("native_compare") !=
                NATIVE_COMPARE_PROTOCOL or
                source_receipt.get("seed_derivation") !=
                SEED_DERIVATION_PROTOCOL or
                source_receipt.get("completion_regime") !=
                provenance["completion_regime"] or
                source_receipt.get("compare_arm") != PRODUCTION_COMPARE_ARM or
                not isinstance(
                    source_receipt.get("measurement_policy"), dict) or
                set(source_receipt["measurement_policy"]) !=
                {"target_profile", "seed_policy", "loss", "schedule"} or
                source_receipt["measurement_policy"].get("target_profile") !=
                TARGET_PROFILE or
                source_receipt["measurement_policy"].get("seed_policy") !=
                TARGET_SEED_POLICY or
                not isinstance(
                    source_receipt["measurement_policy"].get("schedule"),
                    str) or
                source_receipt["measurement_policy"].get("schedule") not in
                TARGET_SCHEDULES or
                type(source_receipt["measurement_policy"].get("loss")) is not
                    float or
                not valid_loss_rate(
                    source_receipt["measurement_policy"].get("loss")) or
                source_receipt.get("recovery_metric_scope") !=
                RECOVERY_METRIC_SCOPE or
                not isinstance(source_receipt.get("python_runtime"), dict) or
                any(
                    type(source_receipt.get(name)) is not int or
                    not 0 <= source_receipt[name] <= 0xffffffffffffffff
                    for name in (
                        "construction_seed_base", "loss_seed_base"))):
            raise MeasurementError(
                "validated peel table has invalid source provenance")
        if (source_benchmark["sha256"] != benchmark["sha256"] or
                source_benchmark["size"] != benchmark["size"]):
            raise MeasurementError(
                "validated peel table changed benchmark identity")
        validation_settings = provenance["settings"]
        required_validation_settings = {
            "source_table", "source_table_sha256", "source_entry_count",
            "selected_entry_count", "selected_K", "trials", "block_bytes",
            "kmax", "margin_percent", "target_profile", "seed_policy",
            "loss", "schedule",
        }
        if (set(validation_settings) != required_validation_settings or
                not isinstance(
                    validation_settings.get("source_table"), str) or
                not validation_settings["source_table"] or
                validation_settings.get("source_table_sha256") !=
                source_provenance["document_sha256"] or
                type(validation_settings.get("source_entry_count")) is not
                    int or
                validation_settings.get("source_entry_count") !=
                source_provenance["entry_count"] or
                type(validation_settings.get("selected_entry_count")) is not
                    int or
                validation_settings.get("selected_entry_count") !=
                source_provenance["selected_entry_count"] or
                not isinstance(
                    validation_settings.get("selected_K"), list) or
                any(type(k) is not int
                    for k in validation_settings["selected_K"]) or
                validation_settings.get("selected_K") !=
                source_provenance["selected_K"] or
                type(validation_settings.get("trials")) is not int or
                not 1 <= validation_settings["trials"] <=
                _COMPARE_TRIALS_MAX or
                type(validation_settings.get("block_bytes")) is not int or
                not 1 <= validation_settings["block_bytes"] <=
                _COMPARE_BLOCK_BYTES_MAX or
                validation_settings["block_bytes"] % 2 != 0 or
                type(validation_settings.get("kmax")) is not int or
                validation_settings["kmax"] < 2 or
                type(validation_settings.get("margin_percent")) is not float or
                not math.isfinite(
                    float(validation_settings["margin_percent"])) or
                validation_settings["margin_percent"] < 0.0 or
                type(validation_settings.get("loss")) is not float or
                validation_settings["loss"] !=
                measurement_policy["loss"] or
                validation_settings.get("target_profile") !=
                TARGET_PROFILE or
                validation_settings.get("seed_policy") !=
                TARGET_SEED_POLICY or
                validation_settings.get("schedule") !=
                measurement_policy["schedule"]):
            raise MeasurementError(
                "validated peel table settings contradict its source receipt")
        expected_search_construction_seed_base = (
            source_receipt["construction_seed_base"])
        expected_search_loss_seed_base = source_receipt["loss_seed_base"]
        source_generator = source_receipt.get("generator")
        source_settings = source_receipt["settings"]
        source_measurement_policy = source_receipt["measurement_policy"]
    else:
        source_measurement_policy = measurement_policy
    search_protocols = {
        "tools/peel_direct.py": (
            "direct-search", "direct-real-codec"),
        "tools/peel_sweep.py": (
            "funnel-search", "unverified-proxy-funnel"),
    }
    if source_generator not in search_protocols:
        raise MeasurementError(
            f"peel table has unsupported search generator "
            f"{source_generator!r}")
    _validate_search_generator_settings(
        source_generator, source_settings, source_measurement_policy,
        "peel table search settings")
    expected_seed_domain, expected_search_kind = search_protocols[
        source_generator]
    if (source_generator == "tools/peel_sweep.py" and
            source_settings.get("allow_unverified_cost_model") is not True):
        raise MeasurementError(
            "proxy funnel table is missing its unverified-cost opt-in")
    raw_entries = document.get("entries")
    if not isinstance(raw_entries, dict) or not raw_entries:
        raise MeasurementError("peel table entries must be a non-empty object")
    entries = {}
    for key, value in raw_entries.items():
        try:
            block_count = int(key)
        except (TypeError, ValueError):
            raise MeasurementError(f"invalid peel table K key: {key!r}")
        if not 2 <= block_count <= 64000 or str(block_count) != key:
            raise MeasurementError(f"non-canonical peel table K key: {key!r}")
        if (not isinstance(value, dict) or
                type(value.get("K")) is not int or
                value.get("K") != block_count):
            raise MeasurementError(
                f"peel table entry {key} has an invalid K receipt")
        coordinate_names = ("scale", "p1", "tilt", "dmax", "absorb")
        if any(name not in value for name in coordinate_names):
            raise MeasurementError(
                f"peel table entry {key} is missing peel coordinates")
        try:
            scale = value["scale"]
        except (KeyError, TypeError):
            raise MeasurementError(
                f"peel table entry {key} has non-numeric coordinates")
        p1 = value["p1"]
        tilt = value["tilt"]
        dmax = value["dmax"]
        absorb = value["absorb"]
        if (type(value["scale"]) is not float or
                any(type(item) is not int
                    for item in (p1, tilt, dmax, absorb)) or
                not math.isfinite(scale) or
                (scale != -1.0 and not 0.0 <= scale <= 64000.0) or
                not 0 <= p1 <= 400 or
                not -100 <= tilt <= 1600 or
                not 2 <= dmax <= 64 or
                not 0 <= absorb <= 100 or
                ("reverted_to_shipped" in value and
                 not isinstance(value["reverted_to_shipped"], bool)) or
                "reverted_to_control" in value):
            raise MeasurementError(
                f"peel table entry {key} has out-of-domain coordinates")
        if value.get("reverted_to_shipped") and (
                scale != -1.0 or p1 != 100 or tilt != 0 or
                dmax != 64 or absorb != 100):
            raise MeasurementError(
                f"peel table entry {key} has a malformed shipped receipt")
        native = value.get("native")
        _validate_native_receipt(
            native, block_count, f"peel table entry {key}.native")
        selected_pmf_digest = _validate_pmf(
            value.get("peel_pmf"),
            f"peel table entry {key}.peel_pmf")
        _validate_search_receipt(
            value.get("search_receipt"), block_count,
            expected_search_construction_seed_base,
            expected_search_loss_seed_base, expected_seed_domain,
            expected_search_kind, native,
            (source_receipt["measurement_policy"]
             if provenance["generator"] == "tools/peel_validate.py"
             else measurement_policy),
            f"peel table entry {key}.search_receipt")
        _validate_entry_generator_settings(
            source_generator, source_settings, value["search_receipt"],
            block_count, f"peel table entry {key}.search_receipt")
        validation_receipt = value.get("validation_receipt")
        if provenance["generator"] == "tools/peel_validate.py":
            _validate_validation_receipt(
                validation_receipt, block_count,
                provenance["construction_seed_base"],
                provenance["loss_seed_base"], native, measurement_policy,
                f"peel table entry {key}.validation_receipt")
            if (validation_receipt["trials"] !=
                    validation_settings["trials"] or
                    validation_receipt["block_bytes"] !=
                    validation_settings["block_bytes"] or
                    validation_receipt["margin_percent"] !=
                    validation_settings["margin_percent"] or
                    validation_receipt["scale"] !=
                    value["search_receipt"]["coordinates"]["scale"] or
                    validation_receipt["trained_pmf_sha256"] != (
                        STOCK_PMF_DIGEST
                        if value["search_receipt"]["peel_pmf"] ==
                        native["pmf"] else
                        pmf_sha256(
                            value["search_receipt"]["peel_pmf"]))):
                raise MeasurementError(
                    f"peel table entry {key} validation receipt contradicts "
                    f"its provenance settings")
            selected_shipped = validation_receipt["verdict"] == "control"
            search_mode = value["search_receipt"]["mode"]
            if (selected_shipped !=
                    bool(value.get("reverted_to_shipped"))):
                raise MeasurementError(
                    f"peel table entry {key} contradicts its validation verdict")
            trained_goodput = _canonical_goodput(
                validation_receipt["trained"], block_count)
            shipped_goodput = _canonical_goodput(
                validation_receipt["shipped"], block_count)
            candidate_clears_margin = (
                search_mode in ("trained", "scale-only") and
                trained_goodput >
                shipped_goodput *
                (1.0 + validation_receipt["margin_percent"] / 100.0)
            )
            if candidate_clears_margin == selected_shipped:
                raise MeasurementError(
                    f"peel table entry {key} validation verdict is inconsistent")
            expected_selected_pmf = (
                native["pmf"] if selected_shipped else
                value["search_receipt"]["peel_pmf"])
            selected_metrics = validation_receipt[
                "shipped" if selected_shipped else "trained"]
            selected_goodput = (
                shipped_goodput if selected_shipped else trained_goodput)
            selected_trials = validation_receipt["trials"]
            validation_delta = (
                None if shipped_goodput == 0.0 else
                round(
                    100.0 * (trained_goodput - shipped_goodput) /
                    shipped_goodput,
                    2)
            )
            expected_gain = validation_delta if not selected_shipped else 0.0
            expected_summaries = {
                "verified_mbps": selected_metrics["decode_mbps"],
                "verified_oh": selected_metrics["oh_mean"],
                "shipped_mbps":
                    validation_receipt["shipped"]["decode_mbps"],
                "gain_pct": expected_gain,
            }
            if selected_shipped:
                expected_summaries[
                    "search_would_have_lost_pct"] = validation_delta
            if (any(
                    name not in value or
                    (expected is None and value[name] is not None) or
                    (expected is not None and (
                        type(value[name]) is not float or
                        value[name] != expected))
                    for name, expected in expected_summaries.items()) or
                    (not selected_shipped and
                     "search_would_have_lost_pct" in value)):
                raise MeasurementError(
                    f"peel table entry {key} has stale validation summary")
        elif validation_receipt is not None:
            raise MeasurementError(
                f"peel table entry {key} has an unexpected validation receipt")
        else:
            search_shipped = value["search_receipt"]["mode"] == "shipped"
            if search_shipped != bool(value.get("reverted_to_shipped")):
                raise MeasurementError(
                    f"peel table entry {key} contradicts its search receipt")
            expected_selected_pmf = value["search_receipt"]["peel_pmf"]
            selected_metrics = value["search_receipt"]
            selected_goodput = _canonical_goodput(
                selected_metrics, block_count)
            selected_trials = value["search_receipt"]["trials"]
        _validate_entry_field_schema(
            value, source_generator, native, value["search_receipt"],
            validated=provenance["generator"] == "tools/peel_validate.py",
            selected_shipped=(
                selected_shipped
                if provenance["generator"] == "tools/peel_validate.py"
                else search_shipped),
            label=f"peel table entry {key}",
        )
        _validate_recovery_metrics(
            value, selected_trials, f"peel table entry {key} top-level")
        if type(value.get("goodput")) is not float:
            raise MeasurementError(
                f"peel table entry {key} has invalid top-level goodput")
        if selected_metrics["fail"] != 0:
            raise MeasurementError(
                f"peel table entry {key} selected a non-decoding arm")
        if (value["peel_pmf"] != expected_selected_pmf or
                selected_pmf_digest != pmf_sha256(expected_selected_pmf)):
            raise MeasurementError(
                f"peel table entry {key} selected PMF is inconsistent")
        top_coordinates = {name: value[name] for name in coordinate_names}
        expected_top_coordinates = value["search_receipt"]["coordinates"]
        if (provenance["generator"] == "tools/peel_validate.py" and
                selected_shipped):
            expected_top_coordinates = {
                "scale": -1.0,
                "p1": 100,
                "tilt": 0,
                "dmax": 64,
                "absorb": 100,
            }
        if top_coordinates != expected_top_coordinates:
            raise MeasurementError(
                f"peel table entry {key} coordinates contradict its search")
        metric_names = (
            "construction_seed", "loss_seed", "target_receipt", "fail",
            "oh_mean", "OH_sd", "OH50", "OH95", "OH99", "OH_max",
            "decode_mbps",
        )
        if (any(value.get(name) != selected_metrics[name]
                for name in metric_names) or
                value.get("goodput") != selected_goodput):
            raise MeasurementError(
                f"peel table entry {key} top-level metrics are stale")
        entries[block_count] = value
    complete_k, expected_selected_k = _expected_source_k(
        source_generator, source_settings, source_selection,
        validation_settings)
    selected_k = sorted(entries)
    if selected_k != expected_selected_k:
        raise MeasurementError(
            "peel table K coverage contradicts its generator settings")
    if source_selection is not None and (
            source_selection["entry_count"] != len(complete_k) or
            source_selection["selected_entry_count"] !=
                len(expected_selected_k) or
            source_selection["selected_K"] != expected_selected_k):
        raise MeasurementError(
            "validated peel table source selection contradicts its entries")
    return entries


def validate_table_document(document):
    """Validate a current table, normalizing numeric overflow to refusal."""
    try:
        return _validate_table_document_impl(document)
    except (OverflowError, ValueError) as error:
        # JSON integer tokens have arbitrary precision in Python.  Every
        # schema numeric eventually enters binary64 or a fixed native integer
        # domain; an integer too large for binary64 conversion raises
        # OverflowError, while one too large for Python's guarded decimal
        # conversion can raise ValueError while deriving a receipt seed.
        # Both are malformed input, not unexpected validator crashes.
        raise MeasurementError(
            f"peel table contains an out-of-range numeric value: {error}")


def _strict_json_object(pairs):
    """Reject duplicate object keys instead of silently accepting the last."""
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value):
    raise ValueError(f"non-standard JSON constant {value!r}")


def _strict_json_float(value):
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite JSON number {value!r}")
    return parsed


def strict_json_loads(payload):
    """Decode strict UTF-8 JSON with unique keys and finite numbers."""
    if isinstance(payload, bytes):
        payload = payload.decode("utf-8")
    if not isinstance(payload, str):
        raise ValueError("JSON payload must be text or bytes")
    return json.loads(
        payload,
        object_pairs_hook=_strict_json_object,
        parse_float=_strict_json_float,
        parse_constant=_reject_json_constant)


def read_table_document_snapshot(path):
    """Read, hash, parse, and validate one stable regular-file snapshot."""
    try:
        with open(path, "rb") as source:
            before = os.fstat(source.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise MeasurementError(
                    f"peel table is not a regular file: {path!r}")
            payload = source.read()
            after = os.fstat(source.fileno())
    except OSError as error:
        raise MeasurementError(f"could not read peel table {path!r}: {error}")
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns")
    if (not stat.S_ISREG(before.st_mode) or
            any(getattr(before, name) != getattr(after, name)
                for name in stable_fields) or
            len(payload) != after.st_size):
        raise MeasurementError(
            f"peel table changed while it was being read: {path!r}")
    try:
        document = strict_json_loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise MeasurementError(f"invalid peel table JSON: {error}")
    entries = validate_table_document(document)
    return document, entries, hashlib.sha256(payload).hexdigest()


def read_table_document(path):
    """Read a current peel table, refusing legacy unversioned results."""
    document, entries, _ = read_table_document_snapshot(path)
    return document, entries


def write_json_atomic(path, value):
    """Replace a JSON result only after the complete value is available."""
    directory = os.path.dirname(os.path.abspath(path))
    handle = tempfile.NamedTemporaryFile(
        mode="w", dir=directory, prefix=".peel-table-",
        suffix=".tmp", delete=False)
    try:
        with handle:
            json.dump(value, handle, indent=1, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(handle.name, path)
    except BaseException:
        try:
            os.unlink(handle.name)
        except FileNotFoundError:
            pass
        raise
