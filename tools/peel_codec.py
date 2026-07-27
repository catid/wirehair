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


PEEL_TABLE_SCHEMA = "wirehair-v2-peel-table/v2"
NATIVE_PMF_PROTOCOL = "wirehair-v2-bench:peelpmf:v1"
NATIVE_COMPARE_PROTOCOL = "wirehair-v2-bench:compare:v2"
SEED_DERIVATION_PROTOCOL = "blake2b64:wh2-peel-v1:domain-separated"
COMPLETION_REGIME_PROTOCOL = "wirehair-v2-production-mixed-completion/v1"
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
PRODUCTION_COMPARE_ARM = {
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
    "seed", "fail", "oh_mean", "OH_sd", "OH50", "OH95", "OH99",
    "OH_max", "decode_mbps",
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


@dataclass(frozen=True)
class StockProfile:
    """Exact native peel distribution and construction metadata for one K."""

    block_count: int
    staircase: int
    source_hits: int
    shipped_mean: float
    pmf: tuple

    def as_dict(self):
        return {
            "block_count": self.block_count,
            "staircase": self.staircase,
            "source_hits": self.source_hits,
            "shipped_mean": self.shipped_mean,
            "pmf_sha256": pmf_sha256(self.pmf),
            "pmf": list(self.pmf),
        }


@dataclass(frozen=True)
class RecoveryMetrics:
    """One independently seeded real-codec recovery measurement."""

    seed: int
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
            "seed": self.seed,
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
        base_seed, seed_domain, coordinates, peel_pmf, shipped_control,
        shipped_goodput, context=None):
    """Create a complete selected-arm receipt for a search result."""
    receipt = metrics.as_dict()
    receipt.update({
        "mode": mode,
        "goodput": goodput,
        "trials": trials,
        "block_bytes": block_bytes,
        "search_kind": search_kind,
        "base_seed": base_seed,
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
    digest = hashlib.blake2b(digest_size=8, person=b"wh2-peel-v1")
    digest.update(base_seed.to_bytes(8, "little"))
    for item in domain:
        encoded = str(item).encode("utf-8")
        digest.update(len(encoded).to_bytes(4, "little"))
        digest.update(encoded)
    return int.from_bytes(digest.digest(), "little")


_STOCK_CACHE = {}


def stock_profile(bench, block_count):
    """Return the exact native PMF and construction metadata for one K."""
    if (isinstance(block_count, bool) or
            not isinstance(block_count, int) or
            not 2 <= block_count <= 64000):
        raise ValueError("native PMF K must be in [2,64000]")
    benchmark = _benchmark_identity(bench)
    key = (benchmark["path"], benchmark["sha256"], block_count)
    cached = _STOCK_CACHE.get(key)
    if cached is not None:
        return cached

    stdout = _run_checked(
        [bench, "peelpmf", "--N", str(block_count)],
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
                if set(values) != {
                        "N", "degrees", "staircase", "source_hits",
                        "shipped_mean"}:
                    raise ValueError("unexpected metadata fields")
                metadata = {
                    "N": int(values["N"]),
                    "degrees": int(values["degrees"]),
                    "staircase": int(values["staircase"]),
                    "source_hits": int(values["source_hits"]),
                    "shipped_mean": float(values["shipped_mean"]),
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
    expected_shipped_mean = (
        block_count *
        min(metadata["source_hits"], metadata["staircase"]) /
        metadata["staircase"]
        if metadata["staircase"] > 0 else math.nan
    )
    if (metadata["N"] != block_count or metadata["degrees"] != 64 or
            metadata["staircase"] < 1 or
            not 1 <= metadata["source_hits"] <= 8 or
            not math.isfinite(metadata["shipped_mean"]) or
            metadata["shipped_mean"] <= 0.0 or
            not math.isclose(
                metadata["shipped_mean"], expected_shipped_mean,
                rel_tol=1e-15, abs_tol=1e-15)):
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
        staircase=metadata["staircase"],
        source_hits=metadata["source_hits"],
        shipped_mean=metadata["shipped_mean"],
        pmf=pmf,
    )
    _STOCK_CACHE[key] = profile
    return profile


def stock_pmf(bench, block_count):
    """Return degrees 1..64 of the exact PMF recovered by the native codec."""
    return list(stock_profile(bench, block_count).pmf)


def pmf_sha256(pmf):
    """Hash the exact binary64 text law passed through the benchmark hook."""
    try:
        encoded = ",".join(
            f"{float(probability):.17g}" for probability in pmf).encode("ascii")
    except (TypeError, ValueError):
        raise ValueError("invalid PMF")
    return hashlib.sha256(encoded).hexdigest()


def family(stock, p1, tilt, dmax, absorb):
    """Apply the native search family's coordinates to a shipped PMF."""
    try:
        coordinates = tuple(float(value)
                            for value in (p1, tilt, dmax, absorb))
    except (TypeError, ValueError):
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
        degree_scale=None, *, seed):
    """Measure one exact codec arm with an explicit independent uint64 seed."""
    if (isinstance(seed, bool) or not isinstance(seed, int) or
            not 0 <= seed <= 0xffffffffffffffff):
        raise ValueError("compare seed must be a uint64")
    if (isinstance(block_count, bool) or
            not isinstance(block_count, int) or
            not 2 <= block_count <= 64000 or
            isinstance(trials, bool) or
            not isinstance(trials, int) or trials < 1 or
            isinstance(block_bytes, bool) or
            not isinstance(block_bytes, int) or block_bytes < 1):
        raise ValueError("invalid compare K, trial count, or block size")
    overrides = {}
    if peel_weights is not None:
        try:
            numeric_weights = [float(value) for value in peel_weights]
        except (TypeError, ValueError):
            raise ValueError("invalid peel weight vector")
        total_weight = sum(numeric_weights)
        if (not 2 <= len(numeric_weights) <= 64 or
                any(not math.isfinite(value) or value < 0.0
                    for value in numeric_weights) or
                not math.isfinite(total_weight) or not total_weight > 0.0):
            raise ValueError("invalid peel weight vector")
        overrides["WIREHAIR_V2_PEEL_DEGREES"] = ",".join(
            f"{value:.17g}" for value in numeric_weights)
    if degree_scale is not None:
        try:
            numeric_scale = float(degree_scale)
        except (TypeError, ValueError):
            raise ValueError("invalid staircase degree scale")
        if (not math.isfinite(numeric_scale) or
                not 0.0 <= numeric_scale <= 64000.0):
            raise ValueError("invalid staircase degree scale")
        overrides["WIREHAIR_V2_STAIRCASE_DEGREE_SCALE"] = (
            f"{numeric_scale:.17g}")
    command = [
        bench, "compare",
        "--nlo", str(block_count),
        "--nhi", str(block_count),
        "--bb-list", str(block_bytes),
        "--trials", str(trials),
        "--loss", "0.10",
        "--seed", str(seed),
        # A fixed cap can silently suppress the requested row when K * bb is
        # larger than the benchmark's default 128 MiB cap.
        "--max-message-mib", "0",
        "--precode",
        "--precode-profile", "mixed",
    ]
    stdout = _run_checked(command, isolated_codec_env(overrides))
    metadata = None
    header = None
    rows = []
    row_codecs = []
    for line in stdout.splitlines():
        if not line.strip():
            continue
        if line.startswith("# compare:"):
            if metadata is not None or header is not None or row_codecs:
                raise MeasurementError(
                    f"misordered or duplicate compare metadata "
                    f"for K={block_count}")
            metadata = {}
            for field in line[len("# compare:"):].split():
                if "=" not in field:
                    raise MeasurementError(
                        f"malformed compare metadata for K={block_count}: "
                        f"{line}")
                name, value = field.split("=", 1)
                if name in metadata:
                    raise MeasurementError(
                        f"duplicate compare metadata field {name!r} "
                        f"for K={block_count}")
                metadata[name] = value
            continue
        if line.startswith("codec"):
            parsed_header = tuple(line.split())
            if (metadata is None or header is not None or row_codecs or
                    parsed_header != COMPARE_COLUMNS):
                raise MeasurementError(
                    f"unexpected, misordered, or duplicate compare header "
                    f"for K={block_count}")
            header = list(parsed_header)
            continue
        fields = line.split()
        codec = fields[0]
        if codec not in ("baseline", "v2", "v2_mixed"):
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
        if codec == "v2_mixed":
            try:
                values = {
                    name: fields[header.index(name)]
                    for name in (
                        "bb", "trials", "N_mean", "fail", "OH_mean", "OH_sd",
                        "OH50", "OH95", "OH99", "OH_max", "decode_MBps")
                }
                metrics = RecoveryMetrics(
                    seed=seed,
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
    if metadata is None:
        raise MeasurementError(
            f"compare returned no metadata for K={block_count}")
    try:
        metadata_seed = int(metadata["seed"], 0)
        schedule_seed = int(metadata["schedule_seed"], 0)
        metadata_loss = float(metadata["loss"])
        dynamic_matches = (
            metadata["N"] == f"[{block_count},{block_count}]" and
            metadata["trials/bb"] == str(trials) and
            metadata["max_message_mib"] == "0" and
            metadata["schedule"] == "iid" and
            metadata["loss_trace"] == "common-id-v2" and
            metadata_seed == seed and schedule_seed == seed and
            math.isfinite(metadata_loss) and metadata_loss == 0.10
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
        *PRODUCTION_COMPARE_ARM,
    }
    if set(metadata) != expected_metadata:
        raise MeasurementError(
            f"compare returned unexpected metadata fields for "
            f"K={block_count}: {sorted(set(metadata) ^ expected_metadata)}")
    completion = {
        name: metadata.get(name)
        for name in PRODUCTION_COMPLETION_REGIME
    }
    if completion != PRODUCTION_COMPLETION_REGIME:
        raise MeasurementError(
            f"compare did not use production mixed completion for "
            f"K={block_count}: {completion}")
    compare_arm = {
        name: metadata.get(name)
        for name in PRODUCTION_COMPARE_ARM
    }
    if compare_arm != PRODUCTION_COMPARE_ARM:
        raise MeasurementError(
            f"compare did not use the production WH2 arm for "
            f"K={block_count}: {compare_arm}")
    if row_codecs != ["baseline", "v2", "v2_mixed"]:
        raise MeasurementError(
            f"compare returned unexpected row sequence {row_codecs!r} "
            f"for K={block_count}")
    if len(rows) != 1:
        raise MeasurementError(
            f"compare returned {len(rows)} v2_mixed rows for K={block_count}")
    return rows[0]


def make_table_document(
        entries, *, generator, bench, base_seed, settings,
        source_provenance=None, artifact_identity=None):
    """Wrap entries in the strict, replayable peel-table schema."""
    if (not isinstance(generator, str) or not generator or
            not isinstance(settings, dict) or
            isinstance(base_seed, bool) or not isinstance(base_seed, int) or
            not 0 <= base_seed <= 0xffffffffffffffff):
        raise ValueError("invalid table generator, settings, or uint64 seed")
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
            "compare_arm": dict(PRODUCTION_COMPARE_ARM),
            "recovery_metric_scope": dict(RECOVERY_METRIC_SCOPE),
            "python_runtime": _python_runtime_identity(),
            "base_seed": base_seed,
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


def _validate_recovery_metrics(receipt, trials, label):
    if (not isinstance(receipt, dict) or
            not RECOVERY_METRIC_FIELDS.issubset(receipt)):
        raise MeasurementError(f"{label} is missing recovery metrics")
    seed = receipt["seed"]
    fail = receipt["fail"]
    if (isinstance(seed, bool) or not isinstance(seed, int) or
            not 0 <= seed <= 0xffffffffffffffff or
            isinstance(fail, bool) or not isinstance(fail, int) or
            not 0 <= fail <= trials):
        raise MeasurementError(f"{label} has an invalid seed or fail count")
    value_names = (
        "oh_mean", "OH_sd", "OH50", "OH95", "OH99", "OH_max",
        "decode_mbps")
    if any(
            isinstance(receipt[name], bool) or
            not isinstance(receipt[name], (int, float))
            for name in value_names):
        raise MeasurementError(f"{label} has non-numeric recovery metrics")
    values = [float(receipt[name]) for name in value_names]
    if (any(not math.isfinite(value) or value < 0.0 for value in values) or
            not values[2] <= values[3] <= values[4] <= values[5] or
            values[0] > values[5]):
        raise MeasurementError(f"{label} has invalid recovery metrics")


def _validate_search_receipt(
        receipt, block_count, expected_base_seed, expected_seed_domain,
        expected_search_kind, native_pmf, label):
    required = {
        "mode", "goodput", "trials", "block_bytes", "search_kind",
        "base_seed", "seed_domain", "coordinates", "peel_pmf",
        "peel_pmf_sha256", "shipped_control", "shipped_goodput", "context",
    }
    if (not isinstance(receipt, dict) or
            set(receipt) != required | RECOVERY_METRIC_FIELDS):
        raise MeasurementError(f"{label} is incomplete")
    trials = receipt["trials"]
    block_bytes = receipt["block_bytes"]
    if (receipt["mode"] not in ("trained", "shipped") or
            isinstance(trials, bool) or not isinstance(trials, int) or
            trials < 1 or isinstance(block_bytes, bool) or
            not isinstance(block_bytes, int) or block_bytes < 1 or
            receipt["search_kind"] != expected_search_kind or
            isinstance(receipt["base_seed"], bool) or
            not isinstance(receipt["base_seed"], int) or
            receipt["base_seed"] != expected_base_seed or
            not 0 <= receipt["base_seed"] <= 0xffffffffffffffff or
            receipt["seed_domain"] != expected_seed_domain or
            not isinstance(receipt["coordinates"], dict) or
            not isinstance(receipt["context"], dict)):
        raise MeasurementError(f"{label} has invalid search metadata")
    expected_sampling_seed = derive_seed(
        expected_base_seed, expected_seed_domain, block_count, "sampling")
    sampling_seed = receipt["context"].get("sampling_seed")
    if (isinstance(sampling_seed, bool) or
            not isinstance(sampling_seed, int) or
            sampling_seed != expected_sampling_seed):
        raise MeasurementError(
            f"{label} has an invalid derived sampling seed")
    coordinate_names = ("scale", "p1", "tilt", "dmax", "absorb")
    coordinates = receipt["coordinates"]
    if set(coordinates) != set(coordinate_names):
        raise MeasurementError(f"{label} has invalid search coordinates")
    scale_value = coordinates["scale"]
    integer_values = [
        coordinates[name] for name in ("p1", "tilt", "dmax", "absorb")
    ]
    if (isinstance(scale_value, bool) or
            not isinstance(scale_value, (int, float)) or
            any(isinstance(value, bool) or not isinstance(value, int)
                for value in integer_values)):
        raise MeasurementError(f"{label} has non-numeric search coordinates")
    scale = float(scale_value)
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
             receipt["mode"] == "trained" and
             not 0.0 <= scale <= 64000.0)):
        raise MeasurementError(f"{label} has invalid search coordinates")
    peel_pmf_digest = _validate_pmf(
        receipt["peel_pmf"], f"{label}.peel_pmf")
    if receipt["peel_pmf_sha256"] != peel_pmf_digest:
        raise MeasurementError(f"{label} has an invalid peel PMF digest")
    if receipt["mode"] == "shipped":
        expected_pmf = list(native_pmf)
        if (scale != -1.0 or p1 != 100 or tilt != 0 or
                dmax != 64 or absorb != 100):
            raise MeasurementError(
                f"{label} has non-production shipped coordinates")
    else:
        expected_pmf = family(native_pmf, p1, tilt, dmax, absorb)
    if expected_pmf is None or receipt["peel_pmf"] != expected_pmf:
        raise MeasurementError(
            f"{label} peel PMF does not match its coordinates")
    # A native-shaped PMF is still a distinct funnel arm when it supplies an
    # active staircase-degree scale.  Direct search has no scale override, so
    # the same PMF there is exactly the shipped arm and must be canonicalized.
    if (receipt["mode"] == "trained" and scale == -1.0 and
            receipt["peel_pmf"] == list(native_pmf)):
        raise MeasurementError(
            f"{label} labels the shipped PMF as a trained arm")
    if (isinstance(receipt["goodput"], bool) or
            not isinstance(receipt["goodput"], (int, float))):
        raise MeasurementError(f"{label} has invalid goodput")
    goodput = float(receipt["goodput"])
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
    expected_seed = derive_seed(
        expected_base_seed, receipt["seed_domain"], block_count, "rank",
        trials, block_bytes)
    if receipt["seed"] != expected_seed:
        raise MeasurementError(f"{label} has an invalid derived rank seed")
    if receipt["shipped_control"]["seed"] != expected_seed:
        raise MeasurementError(
            f"{label} shipped control did not use the paired rank seed")
    expected_goodput = (
        0.0 if receipt["fail"] != 0 else
        receipt["decode_mbps"] * block_count /
        (block_count + receipt["oh_mean"])
    )
    if not math.isclose(
            goodput, expected_goodput, rel_tol=1e-12, abs_tol=1e-12):
        raise MeasurementError(f"{label} has an inconsistent goodput")
    shipped = receipt["shipped_control"]
    expected_shipped_goodput = (
        0.0 if shipped["fail"] != 0 else
        shipped["decode_mbps"] * block_count /
        (block_count + shipped["oh_mean"])
    )
    if (isinstance(receipt["shipped_goodput"], bool) or
            not isinstance(receipt["shipped_goodput"], (int, float)) or
            not math.isclose(
                float(receipt["shipped_goodput"]),
                expected_shipped_goodput,
                rel_tol=1e-12, abs_tol=1e-12)):
        raise MeasurementError(
            f"{label} has an inconsistent shipped-control goodput")
    if (receipt["mode"] == "trained" and
            not goodput > float(receipt["shipped_goodput"])):
        raise MeasurementError(
            f"{label} selected a trained arm that did not beat shipped")
    if receipt["mode"] == "shipped" and any(
            receipt[name] != shipped[name]
            for name in (
                "seed", "fail", "oh_mean", "OH_sd", "OH50", "OH95",
                "OH99", "OH_max", "decode_mbps")):
        raise MeasurementError(
            f"{label} shipped selection contradicts its control metrics")


def _validate_native_receipt(receipt, block_count, label):
    required = {
        "block_count", "staircase", "source_hits", "shipped_mean",
        "pmf_sha256", "pmf",
    }
    if not isinstance(receipt, dict) or set(receipt) != required:
        raise MeasurementError(f"{label} is incomplete")
    if (isinstance(receipt["shipped_mean"], bool) or
            not isinstance(receipt["shipped_mean"], (int, float))):
        raise MeasurementError(f"{label} has invalid shipped mean")
    shipped_mean = float(receipt["shipped_mean"])
    pmf_digest = _validate_pmf(receipt["pmf"], f"{label}.pmf")
    expected_shipped_mean = (
        block_count *
        min(receipt["source_hits"], receipt["staircase"]) /
        receipt["staircase"]
        if (isinstance(receipt["staircase"], int) and
            not isinstance(receipt["staircase"], bool) and
            receipt["staircase"] > 0 and
            isinstance(receipt["source_hits"], int) and
            not isinstance(receipt["source_hits"], bool))
        else math.nan
    )
    if (receipt["block_count"] != block_count or
            isinstance(receipt["block_count"], bool) or
            not isinstance(receipt["block_count"], int) or
            isinstance(receipt["staircase"], bool) or
            not isinstance(receipt["staircase"], int) or
            receipt["staircase"] < 1 or
            isinstance(receipt["source_hits"], bool) or
            not isinstance(receipt["source_hits"], int) or
            not 1 <= receipt["source_hits"] <= 8 or
            not math.isfinite(shipped_mean) or shipped_mean <= 0.0 or
            not math.isclose(
                shipped_mean, expected_shipped_mean,
                rel_tol=1e-15, abs_tol=1e-15) or
            len(receipt["pmf"]) != 64 or
            receipt["pmf_sha256"] != pmf_digest):
        raise MeasurementError(f"{label} has invalid native metadata")


def _validate_pmf(pmf, label):
    if not isinstance(pmf, list) or not 2 <= len(pmf) <= 64:
        raise MeasurementError(f"{label} must contain 2..64 probabilities")
    if any(
            isinstance(value, bool) or not isinstance(value, (int, float))
            for value in pmf):
        raise MeasurementError(f"{label} contains non-numeric probabilities")
    values = [float(value) for value in pmf]
    if (any(not math.isfinite(value) or value < 0.0 for value in values) or
            abs(sum(values) - 1.0) > 1e-12):
        raise MeasurementError(f"{label} is not a probability distribution")
    return pmf_sha256(values)


def _validate_validation_receipt(
        receipt, block_count, base_seed, label):
    required = {
        "verdict", "margin_percent", "trials", "block_bytes",
        "trained_goodput", "shipped_goodput", "trained", "shipped",
    }
    if not isinstance(receipt, dict) or set(receipt) != required:
        raise MeasurementError(f"{label} is incomplete")
    trials = receipt["trials"]
    block_bytes = receipt["block_bytes"]
    if (receipt["verdict"] not in ("keep", "revert") or
            isinstance(trials, bool) or not isinstance(trials, int) or
            trials < 1 or isinstance(block_bytes, bool) or
            not isinstance(block_bytes, int) or block_bytes < 1):
        raise MeasurementError(f"{label} has invalid validation metadata")
    value_names = ("margin_percent", "trained_goodput", "shipped_goodput")
    if any(
            isinstance(receipt[name], bool) or
            not isinstance(receipt[name], (int, float))
            for name in value_names):
        raise MeasurementError(f"{label} has non-numeric validation values")
    values = [float(receipt[name]) for name in value_names]
    if any(not math.isfinite(value) or value < 0.0 for value in values):
        raise MeasurementError(f"{label} has invalid validation values")
    _validate_recovery_metrics(receipt["trained"], trials, f"{label}.trained")
    _validate_recovery_metrics(receipt["shipped"], trials, f"{label}.shipped")
    if (set(receipt["trained"]) != RECOVERY_METRIC_FIELDS or
            set(receipt["shipped"]) != RECOVERY_METRIC_FIELDS):
        raise MeasurementError(f"{label} has unexpected arm fields")
    if receipt["trained"]["seed"] != receipt["shipped"]["seed"]:
        raise MeasurementError(f"{label} did not use a paired compare seed")
    if receipt["trained"]["fail"] > 0 and receipt["shipped"]["fail"] > 0:
        raise MeasurementError(f"{label} records two failed recovery arms")
    expected_seed = derive_seed(
        base_seed, "validation", block_count, trials, block_bytes)
    if receipt["trained"]["seed"] != expected_seed:
        raise MeasurementError(f"{label} has an invalid derived compare seed")
    for arm in ("trained", "shipped"):
        metrics = receipt[arm]
        expected_goodput = (
            0.0 if metrics["fail"] != 0 else
            metrics["decode_mbps"] * block_count /
            (block_count + metrics["oh_mean"])
        )
        if not math.isclose(
                float(receipt[f"{arm}_goodput"]), expected_goodput,
                rel_tol=1e-12, abs_tol=1e-12):
            raise MeasurementError(
                f"{label} has inconsistent {arm} goodput")
    if (receipt["verdict"] == "keep" and
            not receipt["trained_goodput"] >
            receipt["shipped_goodput"] *
            (1.0 + receipt["margin_percent"] / 100.0)):
        raise MeasurementError(
            f"{label} keep verdict does not clear its margin")


def validate_table_document(document):
    """Validate an in-memory current table and return its integer-key entries."""
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
        "seed_derivation", "completion_regime", "base_seed", "settings",
        "source", "compare_arm", "recovery_metric_scope",
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
            isinstance(benchmark.get("size"), bool) or
            not isinstance(benchmark.get("size"), int) or
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
            isinstance(source.get("file_count"), bool) or
            not isinstance(source.get("file_count"), int) or
            source["file_count"] < 1 or
            not isinstance(provenance["settings"], dict) or
            isinstance(provenance["base_seed"], bool) or
            not isinstance(provenance["base_seed"], int) or
            not 0 <= provenance["base_seed"] <= 0xffffffffffffffff):
        raise MeasurementError("peel table has invalid provenance fields")
    expected_search_base_seed = provenance["base_seed"]
    source_generator = provenance["generator"]
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
                isinstance(source_provenance.get("entry_count"), bool) or
                not isinstance(source_provenance.get("entry_count"), int) or
                source_provenance["entry_count"] < 1 or
                isinstance(
                    source_provenance.get("selected_entry_count"), bool) or
                not isinstance(
                    source_provenance.get("selected_entry_count"), int) or
                source_provenance["selected_entry_count"] < 1 or
                not isinstance(source_provenance.get("selected_K"), list) or
                any(isinstance(k, bool) or not isinstance(k, int) or
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
        if (set(source_receipt) != required or
                not isinstance(source_receipt.get("generator"), str) or
                not source_receipt["generator"] or
                not isinstance(source_benchmark, dict) or
                set(source_benchmark) != {"path", "sha256", "size"} or
                not isinstance(source_benchmark.get("path"), str) or
                not source_benchmark["path"] or
                not _is_sha256(source_benchmark.get("sha256")) or
                isinstance(source_benchmark.get("size"), bool) or
                not isinstance(source_benchmark.get("size"), int) or
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
                isinstance(source_identity.get("file_count"), bool) or
                not isinstance(source_identity.get("file_count"), int) or
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
                source_receipt.get("recovery_metric_scope") !=
                RECOVERY_METRIC_SCOPE or
                not isinstance(source_receipt.get("python_runtime"), dict) or
                isinstance(source_receipt.get("base_seed"), bool) or
                not isinstance(source_receipt.get("base_seed"), int) or
                not 0 <= source_receipt["base_seed"] <=
                0xffffffffffffffff):
            raise MeasurementError(
                "validated peel table has invalid source provenance")
        if (source_benchmark["sha256"] != benchmark["sha256"] or
                source_benchmark["size"] != benchmark["size"]):
            raise MeasurementError(
                "validated peel table changed benchmark identity")
        validation_settings = provenance["settings"]
        if (validation_settings.get("source_table_sha256") !=
                source_provenance["document_sha256"] or
                validation_settings.get("source_entry_count") !=
                source_provenance["entry_count"] or
                validation_settings.get("selected_entry_count") !=
                source_provenance["selected_entry_count"] or
                validation_settings.get("selected_K") !=
                source_provenance["selected_K"]):
            raise MeasurementError(
                "validated peel table settings contradict its source receipt")
        expected_search_base_seed = source_receipt["base_seed"]
        source_generator = source_receipt.get("generator")
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
    expected_seed_domain, expected_search_kind = search_protocols[
        source_generator]
    if (provenance["generator"] == "tools/peel_sweep.py" and
            provenance["settings"].get("allow_unverified_cost_model") is
            not True):
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
                isinstance(value.get("K"), bool) or
                not isinstance(value.get("K"), int) or
                value.get("K") != block_count):
            raise MeasurementError(
                f"peel table entry {key} has an invalid K receipt")
        coordinate_names = ("scale", "p1", "tilt", "dmax", "absorb")
        if any(name not in value for name in coordinate_names):
            raise MeasurementError(
                f"peel table entry {key} is missing peel coordinates")
        try:
            scale = float(value["scale"])
        except (OverflowError, TypeError, ValueError):
            raise MeasurementError(
                f"peel table entry {key} has non-numeric coordinates")
        p1 = value["p1"]
        tilt = value["tilt"]
        dmax = value["dmax"]
        absorb = value["absorb"]
        if (isinstance(value["scale"], bool) or
                not isinstance(value["scale"], (int, float)) or
                any(isinstance(item, bool) or not isinstance(item, int)
                    for item in (p1, tilt, dmax, absorb)) or
                not math.isfinite(scale) or
                (scale != -1.0 and not 0.0 <= scale <= 64000.0) or
                not 0 <= p1 <= 400 or
                not -100 <= tilt <= 1600 or
                not 2 <= dmax <= 64 or
                not 0 <= absorb <= 100 or
                ("reverted_to_shipped" in value and
                 not isinstance(value["reverted_to_shipped"], bool))):
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
            expected_search_base_seed, expected_seed_domain,
            expected_search_kind, native["pmf"],
            f"peel table entry {key}.search_receipt")
        validation_receipt = value.get("validation_receipt")
        if provenance["generator"] == "tools/peel_validate.py":
            _validate_validation_receipt(
                validation_receipt, block_count, provenance["base_seed"],
                f"peel table entry {key}.validation_receipt")
            expects_shipped = validation_receipt["verdict"] == "revert"
            if expects_shipped != bool(value.get("reverted_to_shipped")):
                raise MeasurementError(
                    f"peel table entry {key} contradicts its validation verdict")
            trained_clears_margin = (
                value["search_receipt"]["mode"] == "trained" and
                validation_receipt["trained_goodput"] >
                validation_receipt["shipped_goodput"] *
                (1.0 + validation_receipt["margin_percent"] / 100.0)
            )
            if trained_clears_margin == expects_shipped:
                raise MeasurementError(
                    f"peel table entry {key} validation verdict is inconsistent")
            expected_selected_pmf = (
                native["pmf"] if expects_shipped else
                value["search_receipt"]["peel_pmf"])
            selected_metrics = validation_receipt[
                "shipped" if expects_shipped else "trained"]
            selected_goodput = validation_receipt[
                "shipped_goodput" if expects_shipped else "trained_goodput"]
            validation_delta = (
                None if validation_receipt["shipped_goodput"] == 0.0 else
                round(
                    100.0 * (
                        validation_receipt["trained_goodput"] -
                        validation_receipt["shipped_goodput"]) /
                    validation_receipt["shipped_goodput"],
                    2)
            )
            expected_gain = validation_delta if not expects_shipped else 0.0
            if (value.get("verified_mbps") !=
                    selected_metrics["decode_mbps"] or
                    value.get("verified_oh") != selected_metrics["oh_mean"] or
                    value.get("shipped_mbps") !=
                    validation_receipt["shipped"]["decode_mbps"] or
                    value.get("gain_pct") != expected_gain or
                    (expects_shipped and
                     value.get("search_would_have_lost_pct") !=
                     validation_delta)):
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
            selected_goodput = value["search_receipt"]["goodput"]
        if selected_metrics["fail"] != 0:
            raise MeasurementError(
                f"peel table entry {key} selected a non-decoding arm")
        if (value["peel_pmf"] != expected_selected_pmf or
                selected_pmf_digest != pmf_sha256(expected_selected_pmf)):
            raise MeasurementError(
                f"peel table entry {key} selected PMF is inconsistent")
        top_coordinates = {
            name: value[name] for name in coordinate_names
        }
        if (not value.get("reverted_to_shipped") and
                top_coordinates != value["search_receipt"]["coordinates"]):
            raise MeasurementError(
                f"peel table entry {key} coordinates contradict its search")
        metric_names = (
            "seed", "fail", "oh_mean", "OH_sd", "OH50", "OH95", "OH99",
            "OH_max", "decode_mbps",
        )
        if (any(value.get(name) != selected_metrics[name]
                for name in metric_names) or
                value.get("goodput") != selected_goodput):
            raise MeasurementError(
                f"peel table entry {key} top-level metrics are stale")
        entries[block_count] = value
    if source_selection is not None:
        selected_k = sorted(entries)
        if (source_selection["selected_K"] != selected_k or
                source_selection["selected_entry_count"] != len(selected_k) or
                source_selection["entry_count"] <
                source_selection["selected_entry_count"]):
            raise MeasurementError(
                "validated peel table source selection contradicts its entries")
    return entries


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
