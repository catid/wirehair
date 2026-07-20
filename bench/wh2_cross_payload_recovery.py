#!/usr/bin/env python3
"""Sealed two-phase cross-payload recovery evidence for raw WH2 architectures.

This controller deliberately separates cohort discovery from evaluation:

* discovery scans every K with development-only roots at every graph-affecting
  payload width;
* ``seal-cohort`` forms the union of weak K across widths and both raw arms,
  adds only predeclared controls, and freezes the evaluation ledger; and
* evaluation uses disjoint, precommitted roots and common nested OH0/OH1/OH2
  packet prefixes under IID and hard schedules.

The architecture pair is supplied only when ``prepare`` is run.  This lets the
speed campaign select the pair first, while this controller still rejects
K-indexed seed tables, BlockBytes seed normalization, and benchmark-axis flags.
Preparation and reduction are bounded; benchmark work starts only through an
explicit ``run --execute`` command on the frozen controller.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, getcontext
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from typing import (
    Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple,
)


sys.dont_write_bytecode = True
getcontext().prec = 50

SCHEMA = "wirehair.wh2.cross_payload_recovery.v1"
INPUT_SCHEMA = SCHEMA + ".input"
CONTROLLER_NAME = "wh2_cross_payload_recovery.py"
WIDTHS = (64, 256, 1280, 4096)
SCHEDULES = ("iid", "burst", "adversarial", "repair-only")
EVALUATION_LOSSES = ("0.10", "0.35", "0.50")
EVALUATION_OVERHEADS = (0, 1, 2)
DISCOVERY_LOSS = "0.50"
DISCOVERY_OVERHEADS = (0,)
DISCOVERY_TRIALS = 1
EVALUATION_TRIALS = 64
DISCOVERY_ROOT_COUNT = 3
EVALUATION_ROOT_COUNT = 4
DISCOVERY_CHUNK = 256
EVALUATION_CHUNK = 4
K_LO = 2
K_HI = 64000
ARMS = ("control", "candidate")
MAX_STDOUT_BYTES = {
    "discovery": 2 * 1024 * 1024,
    "evaluation": 4 * 1024 * 1024,
}
MAX_STDERR_BYTES = 16 * 1024
# Discovery jobs cover 256 K values x four widths.  Evaluation jobs cover four
# K values x four widths x three nested overheads x 64 trials, exactly three
# times as many trial-solves.  The original 300-second bound was invalidated
# before evaluation after sealed discovery-runtime receipts projected its
# high-K tail past that cap.  1200 seconds is more than 2.5x the projected
# ~470-second tail and changes no seed, cell, outcome, or selection rule.
TIMEOUT_SECONDS = 1200
MAX_WORKERS = 256
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
GIT_OBJECT_RE = re.compile(r"[0-9a-f]{40}\Z")
OPTION_RE = re.compile(r"--[a-z0-9][a-z0-9-]*\Z")

# These fields are printed from wall-clock measurements and are not semantic
# recovery evidence.  The raw stdout remains hash-bound in every receipt.
TIMING_FIELDS = (
    "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
    "residual_ms_mu", "backsub_ms_mu",
)
CSV_FIELDS = (
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

# The controller owns these axes.  Architecture suffixes cannot smuggle a
# second value into a job.  A K-indexed table would make this an architecture
# plus seed-fix comparison, while seed-block normalization would erase the
# width dependence this campaign exists to measure.
OWNED_OPTIONS = {
    "--N", "--bb", "--bb-list", "--overhead", "--trials", "--threads",
    "--loss", "--seed", "--schedule", "--paired-overhead-stream",
    "--full-payload-solve", "--payload-e2e", "--heavy-family",
}
FORBIDDEN_ARCHITECTURE_OPTIONS = {
    "--packet-peel-seed-table", "--seed-block-bytes",
    "--fail-thread-launch-after", "--mixed-null-witnesses",
}
DYNAMIC_PREAMBLE_FIELDS = {
    "trials", "threads", "loss", "seed", "schedule", "overhead_stream",
}
DISCOVERY_SELECTION_INPUTS = (
    "discovery receipts only; evaluation directory absent")
DISCOVERY_SELECTION_RULE = (
    "union of every K with >=1 rank failure in either raw arm at any "
    "width/root/schedule, plus only predeclared controls")


class CampaignError(RuntimeError):
    pass


def die(message: str) -> None:
    raise CampaignError(message)


def canonical_json(value: object) -> bytes:
    return (json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ) + "\n").encode("ascii")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path, maximum: Optional[int] = None) -> str:
    digest = hashlib.sha256()
    total = 0
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1 << 20), b""):
                total += len(block)
                if maximum is not None and total > maximum:
                    die("artifact exceeds its sealed byte bound: {}".format(path))
                digest.update(block)
    except OSError as exc:
        die("cannot hash {}: {}".format(path, exc))
    return digest.hexdigest()


def sealed_record(schema: str, payload: Mapping[str, object]) -> Dict[str, object]:
    if "schema" in payload or "self_sha256_excluding_field" in payload:
        die("sealed payload contains reserved fields")
    record: Dict[str, object] = {"schema": schema, **dict(payload)}
    record["self_sha256_excluding_field"] = sha256_bytes(canonical_json(record))
    return record


def verify_sealed(record: Mapping[str, object], schema: str) -> None:
    if record.get("schema") != schema:
        die("sealed record schema mismatch")
    value = record.get("self_sha256_excluding_field")
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        die("sealed record self hash is malformed")
    payload = dict(record)
    del payload["self_sha256_excluding_field"]
    if sha256_bytes(canonical_json(payload)) != value:
        die("sealed record self hash mismatch")


def write_once(path: Path, data: bytes, mode: int = 0o444) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_symlink():
        die("immutable artifact is a symlink: {}".format(path))
    if path.exists():
        if not path.is_file() or path.read_bytes() != data:
            die("refusing to replace immutable artifact {}".format(path))
        return
    temporary = Path("{}.part.{}.{}".format(
        path, os.getpid(), threading.get_ident()))
    try:
        with temporary.open("xb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(str(temporary), mode)
        try:
            os.link(str(temporary), str(path))
        except FileExistsError:
            if (path.is_symlink() or not path.is_file() or
                    path.read_bytes() != data):
                die("immutable artifact raced with different bytes: {}".format(path))
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def _reject_duplicate_pairs(pairs: Sequence[Tuple[str, object]]) -> Dict[str, object]:
    result: Dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key {}".format(key))
        result[key] = value
    return result


def load_canonical_object(path: Path, context: str) -> Dict[str, object]:
    if not path.is_file() or path.is_symlink():
        die("{} is missing or indirect".format(context))
    try:
        data = path.read_bytes()
        value = json.loads(
            data.decode("ascii"), object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=lambda token: (_ for _ in ()).throw(
                ValueError("nonfinite {}".format(token))),
        )
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as exc:
        die("cannot read {}: {}".format(context, exc))
    if not isinstance(value, dict) or canonical_json(value) != data:
        die("{} is not one canonical JSON object".format(context))
    return value


def load_jsonl(path: Path, context: str) -> List[Dict[str, object]]:
    if not path.is_file() or path.is_symlink():
        die("{} is missing or indirect".format(context))
    try:
        data = path.read_bytes()
    except OSError as exc:
        die("cannot read {}: {}".format(context, exc))
    if not data or not data.endswith(b"\n") or b"\r" in data:
        die("{} is empty or not canonical LF JSONL".format(context))
    result: List[Dict[str, object]] = []
    for index, line in enumerate(data.splitlines(keepends=True)):
        try:
            value = json.loads(
                line.decode("ascii"), object_pairs_hook=_reject_duplicate_pairs,
                parse_constant=lambda token: (_ for _ in ()).throw(
                    ValueError("nonfinite {}".format(token))),
            )
        except (UnicodeError, ValueError, json.JSONDecodeError) as exc:
            die("{} line {} is invalid: {}".format(context, index + 1, exc))
        if not isinstance(value, dict) or canonical_json(value) != line:
            die("{} line {} is not canonical".format(context, index + 1))
        result.append(value)
    return result


def require_plain_int(
    value: object, context: str, minimum: int = 0,
    maximum: int = (1 << 63) - 1,
) -> int:
    if (not isinstance(value, int) or isinstance(value, bool) or
            value < minimum or value > maximum):
        die("{} is outside its integer domain".format(context))
    return value


def require_sha256(value: object, context: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        die("{} is not a lowercase SHA256".format(context))
    return value


def require_git_object(value: object, context: str) -> str:
    if not isinstance(value, str) or GIT_OBJECT_RE.fullmatch(value) is None:
        die("{} is not a full lowercase Git object ID".format(context))
    return value


def require_ascii_text(value: object, context: str, maximum: int = 256) -> str:
    if not isinstance(value, str) or not value or len(value) > maximum:
        die("{} is empty or too long".format(context))
    try:
        value.encode("ascii")
    except UnicodeEncodeError:
        die("{} is not ASCII".format(context))
    if any(ord(char) < 0x20 or ord(char) == 0x7f for char in value):
        die("{} contains control characters".format(context))
    return value


def resolve_direct_path(value: str, context: str) -> Path:
    raw = Path(os.path.abspath(value))
    current = Path(raw.anchor)
    for part in raw.parts[1:]:
        current = current / part
        if current.is_symlink():
            die("{} contains a symlink component: {}".format(context, current))
    return raw.resolve()


def derive_seed(namespace: str, index: int) -> str:
    if index < 0:
        die("negative seed index")
    digest = hashlib.sha256(
        "{}|seed|{}".format(namespace, index).encode("ascii")).digest()
    return "0x{:016x}".format(int.from_bytes(digest[:8], "big"))


DISCOVERY_SEEDS = tuple(
    derive_seed(SCHEMA + ".discovery", index)
    for index in range(DISCOVERY_ROOT_COUNT)
)
EVALUATION_SEEDS = tuple(
    derive_seed(SCHEMA + ".evaluation", index)
    for index in range(EVALUATION_ROOT_COUNT)
)
if set(DISCOVERY_SEEDS).intersection(EVALUATION_SEEDS):
    raise RuntimeError("cross-payload seed namespaces overlap")


def predeclared_control_Ks() -> Tuple[int, ...]:
    """Controls fixed before discovery, including 100 KiB/1 MiB regimes."""
    values: Set[int] = {
        2, 3, 4, 63, 64, 65, 255, 256, 257,
        999, 1000, 1001, 3199, 3200, 3201,
        4095, 4096, 4097, 9999, 10000, 10001,
        19999, 20000, 20001, 31999, 32000, 32001,
        48465, 48466, 48467, 63998, 63999, 64000,
    }
    for width in WIDTHS:
        for target_bytes in (100 * 1024, 1024 * 1024):
            K = (target_bytes + width - 1) // width
            for candidate in (K - 1, K, K + 1):
                if K_LO <= candidate <= K_HI:
                    values.add(candidate)
    return tuple(sorted(values))


CONTROL_KS = predeclared_control_Ks()


def parse_architecture_argv(value: object, context: str) -> Tuple[str, ...]:
    if not isinstance(value, list) or not value:
        die("{} must be a nonempty string list".format(context))
    argv: List[str] = []
    options: Set[str] = set()
    for index, token_value in enumerate(value):
        token = require_ascii_text(token_value, "{} token".format(context), 128)
        if token.startswith("--"):
            if OPTION_RE.fullmatch(token) is None:
                die("{} contains a malformed option {}".format(context, token))
            if token in options:
                die("{} repeats option {}".format(context, token))
            options.add(token)
            if token in OWNED_OPTIONS or token in FORBIDDEN_ARCHITECTURE_OPTIONS:
                die("{} contains controller-owned/forbidden option {}".format(
                    context, token))
            lowered = token.lower()
            if "seed-table" in lowered or "k-fix" in lowered or "fixup" in lowered:
                die("{} contains a K-indexed seed-fix option {}".format(
                    context, token))
        elif index == 0 or not argv[-1].startswith("--"):
            die("{} has a value without a preceding option".format(context))
        argv.append(token)
    try:
        completion = argv[argv.index("--completion") + 1]
        mix_count = argv[argv.index("--mix-count") + 1]
    except (ValueError, IndexError):
        die("{} must explicitly select mixed completion and one mix count".format(
            context))
    if completion != "mixed" or re.fullmatch(r"[1-9][0-9]*", mix_count) is None:
        die("{} does not select one raw mixed-completion architecture".format(
            context))
    return tuple(argv)


def architecture_option_value(argv: Sequence[str], option: str) -> str:
    try:
        value = argv[argv.index(option) + 1]
    except (ValueError, IndexError):
        die("architecture is missing a value for {}".format(option))
    if value.startswith("--"):
        die("architecture is missing a value for {}".format(option))
    return value


def verify_build_receipt_binding(
    receipt: Path, binary: Path, source_commit: str, source_tree: str,
    context: str,
) -> None:
    """Require the copied executable to be a fresh build of the named tree."""
    receipt_record = load_canonical_object(receipt, context)
    source = receipt_record.get("source")
    if (not isinstance(source, dict) or
            source.get("commit") != source_commit or
            source.get("tree_oid") != source_tree or
            receipt_record.get("binary_sha256") != sha256_file(binary) or
            receipt_record.get("fresh_build") is not True):
        die("{} does not bind its source and binary".format(context))


@dataclass(frozen=True)
class ArmInput:
    name: str
    label: str
    binary: Path
    source_commit: str
    source_tree: str
    architecture_argv: Tuple[str, ...]
    build_receipt: Path
    build_receipt_sha256: str

    @staticmethod
    def from_record(record: Mapping[str, object], expected_name: str) -> "ArmInput":
        expected = {
            "name", "label", "binary", "source_commit", "source_tree",
            "architecture_argv", "build_receipt", "build_receipt_sha256",
        }
        if set(record) != expected or record.get("name") != expected_name:
            die("{} arm input fields/name mismatch".format(expected_name))
        binary_value = require_ascii_text(
            record.get("binary"), "{} binary".format(expected_name), 4096)
        receipt_value = require_ascii_text(
            record.get("build_receipt"),
            "{} build receipt".format(expected_name), 4096)
        binary = resolve_direct_path(
            binary_value, "{} binary".format(expected_name))
        receipt = resolve_direct_path(
            receipt_value, "{} build receipt".format(expected_name))
        for path, context in (
                (binary, "binary"), (receipt, "build receipt")):
            if not path.is_file() or path.is_symlink():
                die("{} {} is not a direct regular file".format(
                    expected_name, context))
        if not os.access(str(binary), os.X_OK):
            die("{} binary is not executable".format(expected_name))
        expected_receipt_hash = require_sha256(
            record.get("build_receipt_sha256"),
            "{} build receipt hash".format(expected_name))
        if sha256_file(receipt) != expected_receipt_hash:
            die("{} build receipt hash mismatch".format(expected_name))
        source_commit = require_git_object(
            record.get("source_commit"), "{} source commit".format(expected_name))
        source_tree = require_git_object(
            record.get("source_tree"), "{} source tree".format(expected_name))
        verify_build_receipt_binding(
            receipt, binary, source_commit, source_tree,
            "{} build receipt".format(expected_name))
        return ArmInput(
            name=expected_name,
            label=require_ascii_text(
                record.get("label"), "{} label".format(expected_name), 512),
            binary=binary,
            source_commit=source_commit,
            source_tree=source_tree,
            architecture_argv=parse_architecture_argv(
                record.get("architecture_argv"),
                "{} architecture argv".format(expected_name)),
            build_receipt=receipt,
            build_receipt_sha256=expected_receipt_hash,
        )


def load_input_design(path: Path) -> Tuple[Dict[str, object], Tuple[ArmInput, ...]]:
    record = load_canonical_object(path, "cross-payload input design")
    if set(record) != {"schema", "campaign_label", "arms"} or \
            record.get("schema") != INPUT_SCHEMA:
        die("cross-payload input design fields/schema mismatch")
    require_ascii_text(record.get("campaign_label"), "campaign label", 512)
    arms_value = record.get("arms")
    if not isinstance(arms_value, list) or len(arms_value) != len(ARMS) or \
            any(not isinstance(value, dict) for value in arms_value):
        die("input design must contain exactly two arm objects")
    arms = tuple(
        ArmInput.from_record(arms_value[index], name)
        for index, name in enumerate(ARMS)
    )
    return record, arms


def parse_preamble(line: str) -> Dict[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        die("missing precodefail preamble")
    result: Dict[str, str] = {}
    for token in line[len(prefix):].split():
        if "=" not in token:
            die("malformed precodefail preamble token")
        key, value = token.split("=", 1)
        if not key or not value or key in result:
            die("duplicate/empty precodefail preamble token")
        result[key] = value
    return result


def decimal_field(row: Mapping[str, str], key: str) -> Decimal:
    try:
        value = Decimal(row[key])
    except (KeyError, TypeError, InvalidOperation) as exc:
        die("invalid decimal field {}: {}".format(key, exc))
    if not value.is_finite():
        die("nonfinite decimal field {}".format(key))
    return value


def parse_histogram(text: str, trials: int, context: str) -> Dict[int, int]:
    if not text:
        die("{} histogram is empty".format(context))
    result: Dict[int, int] = {}
    for token in text.split("|"):
        parts = token.split(":")
        if len(parts) != 2 or not all(re.fullmatch(r"[0-9]+", part) for part in parts):
            die("{} histogram is malformed".format(context))
        key, count = int(parts[0]), int(parts[1])
        if count <= 0 or key in result or (result and key <= max(result)):
            die("{} histogram is duplicate/noncanonical".format(context))
        result[key] = count
    if sum(result.values()) != trials:
        die("{} histogram count differs from trials".format(context))
    return result


def parse_failure_trials(text: str, trials: int, count: int) -> Tuple[int, ...]:
    if not text:
        if count:
            die("failure trial list is empty despite failures")
        return ()
    if not all(re.fullmatch(r"[0-9]+", token) for token in text.split("|")):
        die("failure trial list is malformed")
    values = tuple(int(token) for token in text.split("|"))
    if (len(values) != count or values != tuple(sorted(set(values))) or
            any(value < 0 or value >= trials for value in values)):
        die("failure trial list is noncanonical or out of range")
    return values


def validate_row(row: Mapping[str, str], task: Mapping[str, object]) -> Dict[str, object]:
    try:
        K = int(row["N"])
        width = int(row["bb"])
        overhead = int(row["overhead"])
        trials = int(row["trials"])
        success = int(row["success"])
        rank_fail = int(row["rank_fail"])
        error = int(row["error"])
        inact_max = int(row["inact_max"])
        binary_def_max = int(row["binary_def_max"])
        heavy_gain_min = int(row["heavy_gain_min"])
        heavy_shortfall = int(row["heavy_shortfall"])
        seed_attempt = int(row["seed_attempt"])
        first_rank_fail = int(row["first_rank_fail"])
        mix_count = int(row["mix_count"])
    except (KeyError, TypeError, ValueError) as exc:
        die("invalid recovery row integer: {}".format(exc))
    if (K not in task["K"] or width not in task["widths"] or
            overhead not in task["overheads"] or trials != task["trials"]):
        die("recovery row lies outside its exact task grid")
    if min(success, rank_fail, error) < 0 or success + rank_fail + error != trials:
        die("recovery outcome counts do not sum to trials")
    if min(inact_max, binary_def_max, heavy_gain_min, heavy_shortfall, seed_attempt) < 0:
        die("recovery counter is negative")
    failures = parse_failure_trials(
        row["failure_trials"], trials, rank_fail + error)
    # ``failure_trials`` contains both NeedMore and hard-error trials, whereas
    # ``first_rank_fail`` names only the first NeedMore trial.  With no codec
    # errors the first failed trial must therefore be exact.  If errors are
    # present, retain enough structure to receipt the output before the phase
    # reducer rejects the codec error explicitly.
    if ((rank_fail == 0 and first_rank_fail != -1) or
            (rank_fail > 0 and first_rank_fail not in failures) or
            (rank_fail > 0 and error == 0 and first_rank_fail != failures[0])):
        die("first-rank-fail receipt differs from failure trials")
    expected_rate = Decimal(rank_fail + error) / Decimal(trials)
    if abs(decimal_field(row, "fail_rate") - expected_rate) > Decimal("0.0000000051"):
        die("failure rate differs from outcome counts")
    binary_hist = parse_histogram(row["binary_def_hist"], trials, "binary deficit")
    heavy_hist = parse_histogram(row["heavy_gain_hist"], trials, "heavy gain")
    if max(binary_hist) != binary_def_max or min(heavy_hist) != heavy_gain_min:
        die("histogram extrema differ from summary columns")
    decimal_keys = (
        "inact_mu", "binary_def_mu", "heavy_gain_mu", "block_xors_mu",
        "block_muladds_mu", *TIMING_FIELDS,
        "mixed_joint_source_xors_mu", "mixed_joint_marginal_xors_mu",
        "mixed_joint_marginal_copies_mu", "mixed_joint_active_deltas_mu",
        "mixed_joint_scratch_bytes_mu", "mixed_dual_source_columns_mu",
    )
    if any(decimal_field(row, key) < 0 for key in decimal_keys):
        die("recovery metric is negative")
    if decimal_field(row, "inact_mu") > Decimal(inact_max) or \
            decimal_field(row, "binary_def_mu") > Decimal(binary_def_max) or \
            decimal_field(row, "heavy_gain_mu") < Decimal(heavy_gain_min):
        die("recovery mean/extreme receipt is inconsistent")
    if heavy_shortfall > rank_fail or row["heavy_family"] != "periodic":
        die("recovery completion family/shortfall receipt is invalid")
    return {
        "K": K, "width": width, "overhead": overhead, "trials": trials,
        "success": success, "rank_fail": rank_fail, "error": error,
        "failures": failures, "mix_count": mix_count,
        "inact_mu": decimal_field(row, "inact_mu"),
        "block_xors_mu": decimal_field(row, "block_xors_mu"),
        "block_muladds_mu": decimal_field(row, "block_muladds_mu"),
        "seed_attempt": seed_attempt,
    }


def normalized_loss(value: str) -> str:
    # WirehairV2Bench receipts the parsed C++ double with %.17g.  Reproduce
    # that representation exactly rather than comparing a decimal spelling.
    return format(float(value), ".17g")


def expected_dynamic_preamble(task: Mapping[str, object]) -> Dict[str, str]:
    return {
        "trials": str(task["trials"]), "threads": "1",
        "loss": normalized_loss(str(task["loss"])), "seed": str(task["seed"]),
        "schedule": str(task["schedule"]),
        "overhead_stream": (
            "paired" if task["phase"] == "evaluation" else "salted"),
    }


def parse_output(
    data: bytes, task: Mapping[str, object], config: Mapping[str, object],
) -> Tuple[List[Dict[str, object]], str]:
    maximum = require_plain_int(
        task.get("stdout_max_bytes"), "task stdout bound", 1,
        MAX_STDOUT_BYTES[str(task["phase"])] )
    if len(data) > maximum:
        die("benchmark stdout exceeds the sealed task bound")
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as exc:
        die("benchmark stdout is not ASCII: {}".format(exc))
    lines = text.splitlines(keepends=True)
    if (not lines or not all(line.endswith("\n") for line in lines) or
            any("\r" in line or line == "\n" for line in lines)):
        die("benchmark stdout is empty, torn, or contains blank records")
    preamble = parse_preamble(lines[0][:-1])
    dynamic = expected_dynamic_preamble(task)
    if any(preamble.get(key) != value for key, value in dynamic.items()):
        die("benchmark dynamic preamble differs from its sealed task")
    static = {key: value for key, value in preamble.items()
              if key not in DYNAMIC_PREAMBLE_FIELDS}
    arm_config = config["arms"][ARMS.index(str(task["arm"]))]
    if static != arm_config["static_preamble"]:
        die("benchmark architecture preamble differs from its frozen arm")
    if (static.get("seed_block_bytes_override") != "0" or
            static.get("packet_peel_seed_table") != "none" or
            static.get("full_payload_solve") != "0"):
        die("benchmark enabled a seed fix/normalization or payload timing mode")
    reader = csv.DictReader(lines[1:])
    if tuple(reader.fieldnames or ()) != CSV_FIELDS:
        die("benchmark CSV header mismatch")
    rows = list(reader)
    if any(None in row or set(row) != set(CSV_FIELDS) or
           any(value is None for value in row.values()) for row in rows):
        die("benchmark CSV has missing or extra fields")
    expected_keys = [
        (K, width, overhead)
        for width in task["widths"] for K in task["K"]
        for overhead in task["overheads"]
    ]
    parsed: List[Dict[str, object]] = []
    observed: List[Tuple[int, int, int]] = []
    for row in rows:
        value = validate_row(row, task)
        observed.append((
            int(value["K"]), int(value["width"]), int(value["overhead"])))
        parsed.append(value)
    if observed != expected_keys:
        die("benchmark rows do not exactly cover the canonical task grid")
    expected_mix = int(arm_config["mix_count"])
    if any(value["mix_count"] != expected_mix for value in parsed):
        die("CSV mix count differs from the architecture preamble")
    if any(row["active_packet_peel_seed_xor"] !=
           static.get("packet_peel_seed_xor") for row in rows):
        die("CSV active packet seed differs from the raw architecture preamble")
    semantic_rows = []
    for row in rows:
        semantic_rows.append({
            key: row[key] for key in CSV_FIELDS if key not in TIMING_FIELDS})
    return parsed, sha256_bytes(canonical_json(semantic_rows))


def task_argv(task: Mapping[str, object], config: Mapping[str, object]) -> List[str]:
    arm = str(task["arm"])
    arm_config = config["arms"][ARMS.index(arm)]
    argv = [
        str(arm_config["binary"]), "precodefail",
        "--bb-list", ",".join(str(value) for value in task["widths"]),
        "--overhead", ",".join(str(value) for value in task["overheads"]),
        "--trials", str(task["trials"]), "--threads", "1",
        "--loss", str(task["loss"]), "--N",
        ",".join(str(value) for value in task["K"]),
        "--seed", str(task["seed"]), "--schedule", str(task["schedule"]),
        *arm_config["architecture_argv"],
    ]
    if task["phase"] == "evaluation":
        argv.append("--paired-overhead-stream")
    return argv


def _chunks(values: Sequence[int], maximum: int) -> Iterable[Tuple[int, ...]]:
    for offset in range(0, len(values), maximum):
        yield tuple(values[offset:offset + maximum])


def _job_record(
    phase: str, job: int, arm: str, root_index: int, seed: str,
    schedule: str, loss: str, trials: int, overheads: Tuple[int, ...],
    Ks: Tuple[int, ...], cohort_sha256: Optional[str],
    config: Mapping[str, object],
) -> Dict[str, object]:
    record: Dict[str, object] = {
        "schema": SCHEMA + ".job", "phase": phase, "job": job,
        "arm": arm, "root_index": root_index, "seed": seed,
        "schedule": schedule, "loss": loss, "trials": trials,
        "widths": list(WIDTHS), "overheads": list(overheads), "K": list(Ks),
        "cohort_sha256": cohort_sha256,
        "stdout_max_bytes": MAX_STDOUT_BYTES[phase],
        "stderr_max_bytes": MAX_STDERR_BYTES,
        "timeout_seconds": TIMEOUT_SECONDS,
    }
    record["argv"] = task_argv(record, config)
    return record


def build_discovery_jobs(
    config: Mapping[str, object], *, k_lo: int = K_LO, k_hi: int = K_HI,
    chunk: int = DISCOVERY_CHUNK,
) -> List[Dict[str, object]]:
    if not K_LO <= k_lo <= k_hi <= K_HI or not 1 <= chunk <= DISCOVERY_CHUNK:
        die("invalid discovery ledger bounds")
    jobs: List[Dict[str, object]] = []
    Ks = tuple(range(k_lo, k_hi + 1))
    for root_index, seed in enumerate(DISCOVERY_SEEDS):
        for schedule in SCHEDULES:
            for group in _chunks(Ks, chunk):
                for arm in ARMS:
                    jobs.append(_job_record(
                        "discovery", len(jobs), arm, root_index, seed,
                        schedule, DISCOVERY_LOSS, DISCOVERY_TRIALS,
                        DISCOVERY_OVERHEADS, group, None, config))
    return jobs


def build_evaluation_jobs(
    config: Mapping[str, object], cohort: Sequence[int], cohort_sha256: str,
    *, chunk: int = EVALUATION_CHUNK,
) -> List[Dict[str, object]]:
    require_sha256(cohort_sha256, "cohort hash")
    values = tuple(cohort)
    if (not values or
            any(not isinstance(K, int) or isinstance(K, bool) or
                K < K_LO or K > K_HI for K in values) or
            values != tuple(sorted(set(values))) or
            not set(CONTROL_KS).issubset(values) or
            not 1 <= chunk <= EVALUATION_CHUNK):
        die("evaluation cohort/chunk is noncanonical")
    jobs: List[Dict[str, object]] = []
    for root_index, seed in enumerate(EVALUATION_SEEDS):
        for schedule in SCHEDULES:
            for loss in EVALUATION_LOSSES:
                for group in _chunks(values, chunk):
                    for arm in ARMS:
                        jobs.append(_job_record(
                            "evaluation", len(jobs), arm, root_index, seed,
                            schedule, loss, EVALUATION_TRIALS,
                            EVALUATION_OVERHEADS, group, cohort_sha256, config))
    return jobs


def validate_jobs(
    jobs: Sequence[Mapping[str, object]], phase: str,
    config: Mapping[str, object], cohort: Optional[Sequence[int]] = None,
    cohort_sha256: Optional[str] = None,
    *, exact: bool = True,
) -> None:
    if phase not in ("discovery", "evaluation"):
        die("unknown campaign phase")
    for index, task in enumerate(jobs):
        if (task.get("schema") != SCHEMA + ".job" or task.get("phase") != phase or
                task.get("job") != index or task.get("arm") not in ARMS or
                task.get("widths") != list(WIDTHS) or
                task.get("stderr_max_bytes") != MAX_STDERR_BYTES or
                task.get("stdout_max_bytes") != MAX_STDOUT_BYTES[phase] or
                task.get("timeout_seconds") != TIMEOUT_SECONDS):
            die("{} job {} fixed fields changed".format(phase, index))
        if task.get("argv") != task_argv(task, config):
            die("{} job {} command differs from regeneration".format(phase, index))
        if index % 2 == 1:
            left, right = jobs[index - 1], task
            if left.get("arm") != "control" or right.get("arm") != "candidate" or \
                    any(left.get(key) != right.get(key) for key in (
                        "root_index", "seed", "schedule", "loss", "trials",
                        "widths", "overheads", "K", "cohort_sha256")):
                die("adjacent jobs are not an exact raw-architecture pair")
    if len(jobs) % 2:
        die("campaign ledger has an unpaired job")
    if not exact:
        return
    expected = (
        build_discovery_jobs(config) if phase == "discovery" else
        build_evaluation_jobs(config, tuple(cohort or ()), str(cohort_sha256))
    )
    if list(jobs) != expected:
        die("{} ledger does not equal exact regeneration".format(phase))


def manifest_bytes(jobs: Sequence[Mapping[str, object]]) -> bytes:
    return b"".join(canonical_json(dict(job)) for job in jobs)


def discovery_plan_record(
    root: Path, jobs: Sequence[Mapping[str, object]],
) -> Dict[str, object]:
    return sealed_record(SCHEMA + ".discovery_plan", {
        "config_sha256": sha256_file(root / "config.json"),
        "task_manifest_sha256": sha256_file(
            root / "discovery" / "tasks.jsonl"),
        "job_count": len(jobs), "task_pairs": len(jobs) // 2,
        "K_range": [K_LO, K_HI], "K_count": K_HI - K_LO + 1,
        "widths": list(WIDTHS), "roots": len(DISCOVERY_SEEDS),
        "schedules": list(SCHEDULES), "loss": DISCOVERY_LOSS,
        "trials_per_cell": DISCOVERY_TRIALS,
        "selection_has_no_evaluation_input": True,
    })


def evaluation_plan_record(
    cohort_sha256: str, task_manifest_sha256: str,
    jobs: Sequence[Mapping[str, object]], cohort_K_count: int,
) -> Dict[str, object]:
    return sealed_record(SCHEMA + ".evaluation_plan", {
        "cohort_sha256": require_sha256(cohort_sha256, "cohort hash"),
        "task_manifest_sha256": require_sha256(
            task_manifest_sha256, "evaluation task manifest hash"),
        "job_count": len(jobs), "task_pairs": len(jobs) // 2,
        "cohort_K_count": cohort_K_count, "widths": list(WIDTHS),
        "roots": len(EVALUATION_SEEDS), "schedules": list(SCHEDULES),
        "losses": list(EVALUATION_LOSSES),
        "overheads": list(EVALUATION_OVERHEADS),
        "trials_per_cell": EVALUATION_TRIALS,
        "common_nested_packet_prefixes": True,
        "speed_claim_valid": False,
    })


def verify_exact_plan(path: Path, expected: Mapping[str, object], context: str) -> None:
    observed = load_canonical_object(path, context)
    verify_sealed(observed, str(expected["schema"]))
    if observed != expected:
        die("{} differs from exact regeneration".format(context))


def _probe_task(arm: str, config: Mapping[str, object]) -> Dict[str, object]:
    return _job_record(
        "discovery", 0, arm, 0, DISCOVERY_SEEDS[0], SCHEDULES[0],
        DISCOVERY_LOSS, DISCOVERY_TRIALS, DISCOVERY_OVERHEADS, (2,), None,
        config,
    )


def run_captured(argv: Sequence[str], cwd: Path, timeout: int) -> Tuple[bytes, bytes]:
    try:
        result = subprocess.run(
            list(argv), cwd=str(cwd), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=timeout, check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        die("bounded architecture probe failed: {}".format(exc))
    if result.returncode:
        die("bounded architecture probe exited {}: {}".format(
            result.returncode, result.stderr[-4096:].decode("utf-8", "replace")))
    if result.stderr:
        die("bounded architecture probe produced stderr")
    return result.stdout, result.stderr


def extract_static_preamble(data: bytes) -> Dict[str, str]:
    try:
        first = data.decode("ascii").splitlines()[0]
    except (UnicodeError, IndexError) as exc:
        die("cannot parse bounded architecture probe: {}".format(exc))
    preamble = parse_preamble(first)
    static = {key: value for key, value in preamble.items()
              if key not in DYNAMIC_PREAMBLE_FIELDS}
    if (static.get("completion") != "mixed" or
            static.get("seed_block_bytes_override") != "0" or
            static.get("packet_peel_seed_table") != "none" or
            static.get("full_payload_solve") != "0"):
        die("architecture probe is not raw natural-width mixed completion")
    return static


def prepare_config(
    root: Path, input_record: Mapping[str, object], arms: Sequence[ArmInput],
    controller: Path,
) -> Dict[str, object]:
    frozen = root / "frozen"
    frozen.mkdir(parents=True)
    copied_controller = frozen / CONTROLLER_NAME
    shutil.copyfile(str(controller), str(copied_controller))
    os.chmod(str(copied_controller), 0o555)
    arm_records: List[Dict[str, object]] = []
    for arm in arms:
        binary_name = "wirehair_v2_bench.{}".format(arm.name)
        receipt_name = "build_receipt.{}.json".format(arm.name)
        binary = frozen / binary_name
        receipt = frozen / receipt_name
        shutil.copyfile(str(arm.binary), str(binary))
        os.chmod(str(binary), 0o555)
        shutil.copyfile(str(arm.build_receipt), str(receipt))
        os.chmod(str(receipt), 0o444)
        arm_records.append({
            "name": arm.name, "label": arm.label,
            "binary": "frozen/" + binary_name,
            "binary_sha256": sha256_file(binary),
            "source_commit": arm.source_commit, "source_tree": arm.source_tree,
            "architecture_argv": list(arm.architecture_argv),
            "mix_count": int(architecture_option_value(
                arm.architecture_argv, "--mix-count")),
            "build_receipt": "frozen/" + receipt_name,
            "build_receipt_sha256": sha256_file(receipt),
            "static_preamble": {},
        })
    provisional: Dict[str, object] = {
        "schema": SCHEMA + ".config", "campaign_label": input_record["campaign_label"],
        "controller": "frozen/" + CONTROLLER_NAME,
        "controller_sha256": sha256_file(copied_controller),
        "input_design_sha256": sha256_bytes(canonical_json(dict(input_record))),
        "arms": arm_records, "widths": list(WIDTHS),
        "schedules": list(SCHEDULES), "discovery_loss": DISCOVERY_LOSS,
        "discovery_trials": DISCOVERY_TRIALS,
        "evaluation_losses": list(EVALUATION_LOSSES),
        "evaluation_overheads": list(EVALUATION_OVERHEADS),
        "evaluation_trials": EVALUATION_TRIALS,
        "discovery_seeds": list(DISCOVERY_SEEDS),
        "evaluation_seeds": list(EVALUATION_SEEDS),
        "predeclared_control_K": list(CONTROL_KS),
        "seed_fixes_applied": False,
        "natural_block_bytes_graph_seeding": True,
        "recovery_only_no_speed_claim": True,
    }
    # Run one bounded probe per exact copied binary to freeze the complete
    # effective architecture preamble.  The dynamic task axes are removed.
    for arm_record in arm_records:
        task = _probe_task(str(arm_record["name"]), provisional)
        argv = task_argv(task, provisional)
        argv[0] = str(root / str(arm_record["binary"]))
        stdout, _stderr = run_captured(argv, root, TIMEOUT_SECONDS)
        arm_record["static_preamble"] = extract_static_preamble(stdout)
        # The parser now verifies the probe's axes, CSV, and architecture.
        task["argv"] = task_argv(task, provisional)
        parse_output(stdout, task, provisional)
    return sealed_record(SCHEMA + ".config", {
        key: value for key, value in provisional.items() if key != "schema"
    })


def verify_config(
    root: Path, *, require_frozen_controller: bool = True,
) -> Dict[str, object]:
    config = load_canonical_object(root / "config.json", "campaign config")
    verify_sealed(config, SCHEMA + ".config")
    expected_top = {
        "schema", "self_sha256_excluding_field", "campaign_label",
        "controller", "controller_sha256", "input_design_sha256", "arms",
        "widths", "schedules", "discovery_loss", "discovery_trials",
        "evaluation_losses", "evaluation_overheads", "evaluation_trials",
        "discovery_seeds", "evaluation_seeds", "predeclared_control_K",
        "seed_fixes_applied", "natural_block_bytes_graph_seeding",
        "recovery_only_no_speed_claim",
    }
    if set(config) != expected_top or config.get("widths") != list(WIDTHS) or \
            config.get("schedules") != list(SCHEDULES) or \
            config.get("discovery_loss") != DISCOVERY_LOSS or \
            config.get("discovery_trials") != DISCOVERY_TRIALS or \
            config.get("evaluation_losses") != list(EVALUATION_LOSSES) or \
            config.get("evaluation_overheads") != list(EVALUATION_OVERHEADS) or \
            config.get("evaluation_trials") != EVALUATION_TRIALS or \
            config.get("discovery_seeds") != list(DISCOVERY_SEEDS) or \
            config.get("evaluation_seeds") != list(EVALUATION_SEEDS) or \
            config.get("predeclared_control_K") != list(CONTROL_KS) or \
            config.get("seed_fixes_applied") is not False or \
            config.get("natural_block_bytes_graph_seeding") is not True or \
            config.get("recovery_only_no_speed_claim") is not True:
        die("campaign config constants/policy changed")
    require_ascii_text(config.get("campaign_label"), "frozen campaign label", 512)
    require_sha256(config.get("input_design_sha256"), "input design hash")
    require_sha256(config.get("controller_sha256"), "frozen controller hash")
    if config.get("controller") != "frozen/" + CONTROLLER_NAME:
        die("frozen controller path changed")
    arms = config.get("arms")
    if not isinstance(arms, list) or len(arms) != 2:
        die("campaign config arm cardinality changed")
    expected_arm_fields = {
        "name", "label", "binary", "binary_sha256", "source_commit",
        "source_tree", "architecture_argv", "mix_count", "build_receipt",
        "build_receipt_sha256", "static_preamble",
    }
    for index, arm in enumerate(arms):
        if not isinstance(arm, dict) or set(arm) != expected_arm_fields or \
                arm.get("name") != ARMS[index]:
            die("campaign config arm fields/order changed")
        require_ascii_text(arm.get("label"), "frozen arm label", 512)
        parse_architecture_argv(
            arm.get("architecture_argv"), "frozen {} argv".format(ARMS[index]))
        if arm.get("mix_count") != int(architecture_option_value(
                arm["architecture_argv"], "--mix-count")):
            die("frozen arm mix count differs from its architecture argv")
        require_git_object(arm.get("source_commit"), "frozen source commit")
        require_git_object(arm.get("source_tree"), "frozen source tree")
        artifact_paths: Dict[str, Path] = {}
        for path_key, hash_key in (
                ("binary", "binary_sha256"),
                ("build_receipt", "build_receipt_sha256")):
            relative = arm.get(path_key)
            expected_relative = "frozen/{}{}".format(
                "wirehair_v2_bench." if path_key == "binary" else
                "build_receipt.",
                ARMS[index] if path_key == "binary" else ARMS[index] + ".json")
            if relative != expected_relative:
                die("frozen arm path is unsafe")
            path = root / relative
            if not path.is_file() or path.is_symlink() or \
                    sha256_file(path) != require_sha256(
                        arm.get(hash_key), "frozen arm artifact hash"):
                die("frozen arm artifact changed")
            artifact_paths[path_key] = path
        if not os.access(str(artifact_paths["binary"]), os.X_OK):
            die("frozen arm binary is not executable")
        verify_build_receipt_binding(
            artifact_paths["build_receipt"], artifact_paths["binary"],
            str(arm["source_commit"]), str(arm["source_tree"]),
            "frozen {} build receipt".format(ARMS[index]))
        static = arm.get("static_preamble")
        if not isinstance(static, dict) or any(
                not isinstance(key, str) or not isinstance(value, str)
                for key, value in static.items()) or \
                any(key in DYNAMIC_PREAMBLE_FIELDS for key in static) or \
                static.get("completion") != "mixed" or \
                static.get("seed_block_bytes_override") != "0" or \
                static.get("packet_peel_seed_table") != "none" or \
                static.get("full_payload_solve") != "0":
            die("frozen architecture preamble is malformed or non-raw")
    controller = root / str(config["controller"])
    if not controller.is_file() or controller.is_symlink() or \
            sha256_file(controller) != config.get("controller_sha256"):
        die("frozen controller changed")
    if require_frozen_controller and Path(__file__).resolve() != controller.resolve():
        die("campaign operations must use the exact frozen controller")
    return config


def phase_paths(root: Path, phase: str, job: int) -> Dict[str, Path]:
    stem = "job{:05d}".format(job)
    base = root / phase
    return {
        "stdout": base / "raw" / (stem + ".csv"),
        "stderr": base / "stderr" / (stem + ".txt"),
        "receipt": base / "receipts" / (stem + ".json"),
    }


def verify_phase_layout(root: Path, phase: str) -> None:
    if phase not in ("discovery", "evaluation"):
        die("unknown phase layout")
    base = root / phase
    for path in (base, base / "raw", base / "stderr", base / "receipts"):
        if not path.is_dir() or path.is_symlink():
            die("{} phase directory is missing or indirect: {}".format(
                phase, path))


def verify_complete_phase_inventory(
    root: Path, phase: str, jobs: Sequence[Mapping[str, object]],
) -> None:
    verify_phase_layout(root, phase)
    allowed_top = {"raw", "stderr", "receipts", "tasks.jsonl", "plan.json"}
    result_path = root / phase / "result.json"
    if phase == "evaluation" and result_path.exists():
        if not result_path.is_file() or result_path.is_symlink():
            die("evaluation result is indirect or not a regular file")
        allowed_top.add("result.json")
    top_entries = list((root / phase).iterdir())
    if {path.name for path in top_entries} != allowed_top:
        die("{} top-level inventory contains missing or unexpected entries".format(
            phase))
    expected = {"job{:05d}".format(int(task["job"])) for task in jobs}
    directories = {
        "raw": {stem + ".csv" for stem in expected},
        "stderr": {stem + ".txt" for stem in expected},
        "receipts": {stem + ".json" for stem in expected},
    }
    for name, names in directories.items():
        entries = list((root / phase / name).iterdir())
        if ({path.name for path in entries} != names or
                any(not path.is_file() or path.is_symlink() for path in entries)):
            die("{} {} inventory is incomplete or contains extras".format(
                phase, name))


def receipt_record(
    task: Mapping[str, object], stdout: bytes, stderr: bytes,
    semantic_sha256: str, started_ns: int, ended_ns: int,
) -> Dict[str, object]:
    return sealed_record(SCHEMA + ".job_receipt", {
        "phase": task["phase"], "job": task["job"], "arm": task["arm"],
        "task_sha256": sha256_bytes(canonical_json(dict(task))),
        "argv": task["argv"], "returncode": 0,
        "started_monotonic_ns": started_ns, "ended_monotonic_ns": ended_ns,
        "elapsed_ns": ended_ns - started_ns,
        "stdout_sha256": sha256_bytes(stdout), "stdout_bytes": len(stdout),
        "stderr_sha256": sha256_bytes(stderr), "stderr_bytes": len(stderr),
        "semantic_sha256": semantic_sha256,
        "timing_fields_are_not_recovery_evidence": list(TIMING_FIELDS),
    })


def verify_receipt(
    root: Path, task: Mapping[str, object], config: Mapping[str, object],
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    paths = phase_paths(root, str(task["phase"]), int(task["job"]))
    for path in paths.values():
        if not path.is_file() or path.is_symlink():
            die("job {} artifact is missing or indirect".format(task["job"]))
    receipt = load_canonical_object(paths["receipt"], "job receipt")
    verify_sealed(receipt, SCHEMA + ".job_receipt")
    expected_fields = {
        "schema", "self_sha256_excluding_field", "phase", "job", "arm",
        "task_sha256", "argv", "returncode", "started_monotonic_ns",
        "ended_monotonic_ns", "elapsed_ns", "stdout_sha256", "stdout_bytes",
        "stderr_sha256", "stderr_bytes", "semantic_sha256",
        "timing_fields_are_not_recovery_evidence",
    }
    if set(receipt) != expected_fields:
        die("job receipt fields changed")
    stdout = paths["stdout"].read_bytes()
    stderr = paths["stderr"].read_bytes()
    parsed, semantic = parse_output(stdout, task, config)
    started = require_plain_int(
        receipt.get("started_monotonic_ns"), "receipt start", 1)
    ended = require_plain_int(receipt.get("ended_monotonic_ns"), "receipt end", 1)
    if (receipt.get("phase") != task["phase"] or
            receipt.get("job") != task["job"] or
            receipt.get("arm") != task["arm"] or
            receipt.get("task_sha256") != sha256_bytes(canonical_json(dict(task))) or
            receipt.get("argv") != task["argv"] or receipt.get("returncode") != 0 or
            ended < started or receipt.get("elapsed_ns") != ended - started or
            receipt.get("stdout_sha256") != sha256_bytes(stdout) or
            receipt.get("stdout_bytes") != len(stdout) or
            receipt.get("stderr_sha256") != sha256_bytes(stderr) or
            receipt.get("stderr_bytes") != len(stderr) or stderr or
            receipt.get("semantic_sha256") != semantic or
            receipt.get("timing_fields_are_not_recovery_evidence") !=
                list(TIMING_FIELDS)):
        die("job receipt/output binding mismatch")
    return receipt, parsed


def verify_cohort_record(
    root: Path, config: Mapping[str, object], *, verify_discovery: bool = False,
) -> Dict[str, object]:
    cohort = load_canonical_object(root / "cohort.json", "sealed cohort")
    verify_sealed(cohort, SCHEMA + ".cohort")
    expected_fields = {
        "schema", "self_sha256_excluding_field", "config_sha256",
        "discovery_manifest_sha256", "discovery_receipt_ledger_sha256",
        "discovery_seed_namespace", "evaluation_seed_namespace",
        "discovery_seeds", "evaluation_seed_commitment_sha256",
        "selection_inputs", "selection_rule", "weak_K_by_width",
        "weak_K_by_arm_width", "failure_cells_by_arm_width",
        "weak_K_multiplicity", "weak_union", "controls", "cohort",
    }
    if (set(cohort) != expected_fields or
            cohort.get("config_sha256") != sha256_file(root / "config.json") or
            cohort.get("discovery_manifest_sha256") != sha256_file(
                root / "discovery" / "tasks.jsonl") or
            cohort.get("discovery_seed_namespace") != SCHEMA + ".discovery" or
            cohort.get("evaluation_seed_namespace") != SCHEMA + ".evaluation" or
            cohort.get("discovery_seeds") != list(DISCOVERY_SEEDS) or
            cohort.get("evaluation_seed_commitment_sha256") != sha256_bytes(
                canonical_json(list(EVALUATION_SEEDS))) or
            cohort.get("selection_inputs") != DISCOVERY_SELECTION_INPUTS or
            cohort.get("selection_rule") != DISCOVERY_SELECTION_RULE):
        die("sealed cohort provenance/policy changed")
    require_sha256(
        cohort.get("discovery_receipt_ledger_sha256"),
        "discovery receipt ledger hash")
    weak_union = cohort.get("weak_union")
    controls = cohort.get("controls")
    values = cohort.get("cohort")
    def canonical_K_list(sequence: object) -> bool:
        return (
            isinstance(sequence, list) and
            all(isinstance(K, int) and not isinstance(K, bool) and
                K_LO <= K <= K_HI for K in sequence) and
            sequence == sorted(set(sequence)))

    for name, sequence in (
            ("weak union", weak_union), ("controls", controls),
            ("cohort", values)):
        if not canonical_K_list(sequence):
            die("sealed {} K list is noncanonical".format(name))
    if (not isinstance(weak_union, list) or not isinstance(controls, list) or
            not isinstance(values, list)):
        die("sealed cohort K lists are malformed")
    if controls != list(CONTROL_KS) or values != sorted(set(weak_union).union(controls)):
        die("sealed cohort was not formed from weak union plus declared controls")
    by_width = cohort.get("weak_K_by_width")
    by_arm_width = cohort.get("weak_K_by_arm_width")
    failure_counts = cohort.get("failure_cells_by_arm_width")
    multiplicity = cohort.get("weak_K_multiplicity")
    width_keys = {str(width) for width in WIDTHS}
    if (not isinstance(by_width, dict) or set(by_width) != width_keys or
            not isinstance(by_arm_width, dict) or set(by_arm_width) != set(ARMS) or
            not isinstance(failure_counts, dict) or
                set(failure_counts) != set(ARMS) or
            not isinstance(multiplicity, dict) or set(multiplicity) != set(ARMS)):
        die("sealed cohort width/arm maps are malformed")
    observed_union: Set[int] = set()
    for width in WIDTHS:
        width_key = str(width)
        width_values = by_width[width_key]
        if (not canonical_K_list(width_values) or
                not set(width_values).issubset(weak_union)):
            die("sealed per-width weak-K list is malformed")
        observed_union.update(width_values)
    if sorted(observed_union) != weak_union:
        die("sealed weak union differs from the per-width union")
    arm_width_sets: Dict[str, Dict[str, Set[int]]] = {
        arm: {} for arm in ARMS}
    for arm in ARMS:
        if (not isinstance(by_arm_width[arm], dict) or
                set(by_arm_width[arm]) != width_keys or
                not isinstance(failure_counts[arm], dict) or
                set(failure_counts[arm]) != width_keys or
                not isinstance(multiplicity[arm], dict)):
            die("sealed per-arm cohort maps are malformed")
        arm_union: Set[int] = set()
        for width in WIDTHS:
            key = str(width)
            arm_values = by_arm_width[arm][key]
            if (not canonical_K_list(arm_values) or
                    not set(arm_values).issubset(by_width[key]) or
                    require_plain_int(
                        failure_counts[arm][key], "failure cell count") <
                        len(arm_values) or
                    int(failure_counts[arm][key]) >
                        len(arm_values) * DISCOVERY_ROOT_COUNT * len(SCHEDULES)):
                die("sealed per-arm/width weak evidence is malformed")
            arm_width_sets[arm][key] = set(arm_values)
            arm_union.update(arm_values)
        for K_text, count in multiplicity[arm].items():
            if (not isinstance(K_text, str) or
                    re.fullmatch(r"[1-9][0-9]*", K_text) is None or
                    int(K_text) not in arm_union or
                    require_plain_int(count, "weak-K multiplicity", 1) >
                        DISCOVERY_ROOT_COUNT * len(SCHEDULES) * len(WIDTHS)):
                die("sealed weak-K multiplicity map is malformed")
        if {int(K) for K in multiplicity[arm]} != arm_union:
            die("sealed weak-K multiplicity keys differ from per-width maps")
        if sum(int(value) for value in failure_counts[arm].values()) != \
                sum(int(value) for value in multiplicity[arm].values()):
            die("sealed failure-cell totals differ from weak-K multiplicity")
    for width in WIDTHS:
        key = str(width)
        if set(by_width[key]) != set().union(
                *(arm_width_sets[arm][key] for arm in ARMS)):
            die("sealed per-width weak K is not the exact arm union")
    if verify_discovery:
        jobs = load_jsonl(
            root / "discovery" / "tasks.jsonl", "discovery task ledger")
        validate_jobs(jobs, "discovery", config)
        records, receipt_ledger = collect_phase_records(
            root, "discovery", config, jobs)
        projected: Dict[
            str, Dict[Tuple[int, int, str, int], Mapping[str, object]]
        ] = {arm: {} for arm in ARMS}
        for arm in ARMS:
            for (K, root_index, schedule, loss, width, overhead), row in \
                    records[arm].items():
                if loss != DISCOVERY_LOSS or overhead != 0:
                    die("discovery row contains an evaluation-only axis")
                projected[arm][(K, root_index, schedule, width)] = row
        derived = derive_cohort(projected)
        for key, value in derived.items():
            if cohort.get(key) != value:
                die("sealed cohort differs from discovery receipt reduction")
        if cohort.get("discovery_receipt_ledger_sha256") != sha256_bytes(
                canonical_json(receipt_ledger)):
            die("sealed cohort receipt ledger hash differs from discovery")
    return cohort


def load_phase_jobs(
    root: Path, phase: str, config: Mapping[str, object],
) -> List[Dict[str, object]]:
    verify_phase_layout(root, phase)
    jobs = load_jsonl(root / phase / "tasks.jsonl", phase + " task ledger")
    if phase == "discovery":
        validate_jobs(jobs, phase, config)
        verify_exact_plan(
            root / "discovery" / "plan.json",
            discovery_plan_record(root, jobs), "discovery plan")
    else:
        cohort = verify_cohort_record(root, config)
        validate_jobs(
            jobs, phase, config, cohort["cohort"],
            sha256_bytes(canonical_json(cohort)))
        verify_exact_plan(
            root / "evaluation" / "plan.json",
            evaluation_plan_record(
                sha256_file(root / "cohort.json"),
                sha256_file(root / "evaluation" / "tasks.jsonl"),
                jobs, len(cohort["cohort"])),
            "evaluation plan")
    return jobs


def _execute_job(
    root: Path, task: Mapping[str, object], config: Mapping[str, object],
) -> Tuple[int, bool]:
    paths = phase_paths(root, str(task["phase"]), int(task["job"]))
    if paths["receipt"].exists():
        verify_receipt(root, task, config)
        return int(task["job"]), True
    if any(path.exists() or path.is_symlink() for path in paths.values()):
        die("unreceipted job artifacts require explicit operator inspection")
    argv = list(task["argv"])
    argv[0] = str(root / argv[0])
    started = time.monotonic_ns()
    process: Optional[subprocess.Popen] = None
    try:
        process = subprocess.Popen(
            argv, cwd=str(root), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True,
        )
        try:
            stdout, stderr = process.communicate(timeout=int(task["timeout_seconds"]))
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            die("job {} timed out".format(task["job"]))
    except OSError as exc:
        die("job {} failed to launch: {}".format(task["job"], exc))
    except BaseException:
        if process is not None and process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        raise
    ended = time.monotonic_ns()
    if process.returncode:
        die("job {} exited {}: {}".format(
            task["job"], process.returncode,
            stderr[-4096:].decode("utf-8", "replace")))
    if len(stderr) > int(task["stderr_max_bytes"]) or stderr:
        die("job {} produced stderr".format(task["job"]))
    parsed, semantic = parse_output(stdout, task, config)
    del parsed
    write_once(paths["stdout"], stdout)
    write_once(paths["stderr"], stderr)
    receipt = receipt_record(task, stdout, stderr, semantic, started, ended)
    write_once(paths["receipt"], canonical_json(receipt))
    verify_receipt(root, task, config)
    return int(task["job"]), False


def run_phase(
    root: Path, phase: str, workers: int, execute: bool,
) -> Dict[str, object]:
    if not 1 <= workers <= MAX_WORKERS:
        die("worker count is outside the supported range")
    config = verify_config(root)
    jobs = load_phase_jobs(root, phase, config)
    complete = 0
    for task in jobs:
        path = phase_paths(root, phase, int(task["job"]))["receipt"]
        if path.exists():
            verify_receipt(root, task, config)
            complete += 1
        elif any(candidate.exists() or candidate.is_symlink()
                 for candidate in phase_paths(
                     root, phase, int(task["job"])).values()):
            die("unreceipted job artifacts require explicit operator inspection")
    if not execute:
        return {"phase": phase, "jobs": len(jobs), "complete": complete,
                "remaining": len(jobs) - complete, "executed": False}
    failures: List[str] = []
    resumed = launched = 0
    remaining = [
        task for task in jobs if not phase_paths(
            root, phase, int(task["job"]))["receipt"].exists()
    ]
    iterator = iter(remaining)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        future_map: Dict[concurrent.futures.Future, int] = {}

        def submit_one() -> bool:
            try:
                task = next(iterator)
            except StopIteration:
                return False
            future_map[pool.submit(_execute_job, root, task, config)] = int(task["job"])
            return True

        for _ in range(min(workers, len(remaining))):
            submit_one()
        abort = False
        while future_map:
            done, _pending = concurrent.futures.wait(
                tuple(future_map), return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                job = future_map.pop(future)
                try:
                    _job, was_resumed = future.result()
                    resumed += int(was_resumed)
                    launched += int(not was_resumed)
                except Exception as exc:  # preserve every completed receipt
                    failures.append("job {}: {}".format(job, exc))
                    abort = True
            if abort:
                for future in future_map:
                    future.cancel()
            else:
                for _ in range(len(done)):
                    if not submit_one():
                        break
    if failures:
        die("phase execution failed: {}".format("; ".join(failures[:8])))
    return {"phase": phase, "jobs": len(jobs), "complete": len(jobs),
            "launched": launched, "resumed": resumed, "executed": True}


def derive_cohort(
    records: Mapping[str, Mapping[Tuple[int, int, str, int], Mapping[str, object]]],
    *, k_lo: int = K_LO, k_hi: int = K_HI,
) -> Dict[str, object]:
    expected = {
        (K, root_index, schedule, width)
        for K in range(k_lo, k_hi + 1)
        for root_index in range(DISCOVERY_ROOT_COUNT)
        for schedule in SCHEDULES for width in WIDTHS
    }
    if set(records) != set(ARMS) or any(set(records[arm]) != expected for arm in ARMS):
        die("discovery cube does not exactly cover K x roots x schedules x widths")
    weak_by_width: Dict[int, Set[int]] = {width: set() for width in WIDTHS}
    weak_by_arm_width: Dict[str, Dict[int, Set[int]]] = {
        arm: {width: set() for width in WIDTHS} for arm in ARMS}
    failure_cells_by_arm_width: Dict[str, Dict[int, int]] = {
        arm: {width: 0 for width in WIDTHS} for arm in ARMS}
    multiplicity: Dict[str, Dict[int, int]] = {arm: {} for arm in ARMS}
    for arm in ARMS:
        for (K, _root, _schedule, width), row in records[arm].items():
            if int(row["error"]):
                die("discovery contains a codec error; cohort is invalid")
            if int(row["rank_fail"]):
                weak_by_width[width].add(K)
                weak_by_arm_width[arm][width].add(K)
                failure_cells_by_arm_width[arm][width] += 1
                multiplicity[arm][K] = multiplicity[arm].get(K, 0) + 1
    weak_union: Set[int] = set().union(*weak_by_width.values())
    controls = {K for K in CONTROL_KS if k_lo <= K <= k_hi}
    cohort = tuple(sorted(weak_union.union(controls)))
    return {
        "weak_K_by_width": {
            str(width): sorted(weak_by_width[width]) for width in WIDTHS},
        "weak_K_by_arm_width": {
            arm: {str(width): sorted(weak_by_arm_width[arm][width])
                  for width in WIDTHS} for arm in ARMS},
        "failure_cells_by_arm_width": {
            arm: {str(width): failure_cells_by_arm_width[arm][width]
                  for width in WIDTHS} for arm in ARMS},
        "weak_K_multiplicity": {
            arm: {str(K): multiplicity[arm][K] for K in sorted(multiplicity[arm])}
            for arm in ARMS},
        "weak_union": sorted(weak_union), "controls": sorted(controls),
        "cohort": list(cohort),
    }


def collect_phase_records(
    root: Path, phase: str, config: Mapping[str, object],
    jobs: Sequence[Mapping[str, object]],
) -> Tuple[
        Dict[str, Dict[Tuple[int, int, str, str, int, int], Dict[str, object]]],
        List[Dict[str, object]],
    ]:
    verify_complete_phase_inventory(root, phase, jobs)
    records: Dict[
        str, Dict[Tuple[int, int, str, str, int, int], Dict[str, object]]
    ] = {arm: {} for arm in ARMS}
    receipt_ledger: List[Dict[str, object]] = []
    for task in jobs:
        receipt, rows = verify_receipt(root, task, config)
        receipt_path = phase_paths(root, phase, int(task["job"]))["receipt"]
        receipt_ledger.append({
            "job": task["job"], "receipt_sha256": sha256_file(receipt_path),
            "semantic_sha256": receipt["semantic_sha256"],
        })
        for row in rows:
            key = (
                int(row["K"]), int(task["root_index"]), str(task["schedule"]),
                str(task["loss"]), int(row["width"]), int(row["overhead"]),
            )
            arm = str(task["arm"])
            if key in records[arm]:
                die("duplicate {} semantic cell {}".format(phase, key))
            records[arm][key] = row
    return records, receipt_ledger


def seal_discovery_cohort(root: Path) -> Dict[str, object]:
    config = verify_config(root)
    jobs = load_phase_jobs(root, "discovery", config)
    evaluation = root / "evaluation"
    cohort_path = root / "cohort.json"
    if evaluation.is_symlink() or cohort_path.is_symlink():
        die("cohort/evaluation publication path is indirect")
    if evaluation.exists() and not cohort_path.is_file():
        die("evaluation was published without its cohort")
    records, receipt_ledger = collect_phase_records(
        root, "discovery", config, jobs)
    projected: Dict[str, Dict[Tuple[int, int, str, int], Mapping[str, object]]] = {
        arm: {} for arm in ARMS}
    for arm in ARMS:
        for (K, root_index, schedule, loss, width, overhead), row in records[arm].items():
            if loss != DISCOVERY_LOSS or overhead != 0:
                die("discovery row contains an evaluation-only axis")
            projected[arm][(K, root_index, schedule, width)] = row
    derived = derive_cohort(projected)
    discovery_manifest = root / "discovery" / "tasks.jsonl"
    receipt_ledger_sha256 = sha256_bytes(canonical_json(receipt_ledger))
    cohort = sealed_record(SCHEMA + ".cohort", {
        "config_sha256": sha256_file(root / "config.json"),
        "discovery_manifest_sha256": sha256_file(discovery_manifest),
        "discovery_receipt_ledger_sha256": receipt_ledger_sha256,
        "discovery_seed_namespace": SCHEMA + ".discovery",
        "evaluation_seed_namespace": SCHEMA + ".evaluation",
        "discovery_seeds": list(DISCOVERY_SEEDS),
        "evaluation_seed_commitment_sha256": sha256_bytes(
            canonical_json(list(EVALUATION_SEEDS))),
        "selection_inputs": DISCOVERY_SELECTION_INPUTS,
        "selection_rule": DISCOVERY_SELECTION_RULE,
        **derived,
    })
    cohort_data = canonical_json(cohort)
    cohort_file_sha = sha256_bytes(cohort_data)
    if cohort_path.exists() and (
            not cohort_path.is_file() or cohort_path.read_bytes() != cohort_data):
        die("published cohort differs from exact discovery reduction")
    eval_jobs = build_evaluation_jobs(config, cohort["cohort"], cohort_file_sha)
    eval_data = manifest_bytes(eval_jobs)
    if evaluation.exists():
        tasks_path = evaluation / "tasks.jsonl"
        if (not evaluation.is_dir() or not tasks_path.is_file() or
                tasks_path.is_symlink() or tasks_path.read_bytes() != eval_data):
            die("published evaluation ledger differs from exact regeneration")
        verify_phase_layout(root, "evaluation")
        expected_plan = evaluation_plan_record(
            cohort_file_sha, sha256_bytes(eval_data), eval_jobs,
            len(cohort["cohort"]))
        verify_exact_plan(
            evaluation / "plan.json", expected_plan, "evaluation plan")
        return cohort
    staging = root / ".evaluation.prepare"
    if staging.is_symlink():
        die("evaluation staging path is indirect")
    if staging.exists():
        if not staging.is_dir():
            die("evaluation staging path is not a directory")
        # No command can execute from this unpublished directory.  It contains
        # deterministic plan bytes only, so a torn prior publication is safe
        # to reconstruct from the just-reduced discovery receipts.
        shutil.rmtree(str(staging))
    for directory in ("raw", "stderr", "receipts"):
        (staging / directory).mkdir(parents=True)
    write_once(staging / "tasks.jsonl", eval_data)
    audit = evaluation_plan_record(
        cohort_file_sha, sha256_file(staging / "tasks.jsonl"), eval_jobs,
        len(cohort["cohort"]))
    write_once(staging / "plan.json", canonical_json(audit))
    write_once(cohort_path, cohort_data)
    os.rename(str(staging), str(evaluation))
    return cohort


def _axis_key(name: str, value: object) -> str:
    return "{}={}".format(name, value)


def summarize_evaluation_records(
    records: Mapping[
        str, Mapping[
            Tuple[int, int, str, str, int, int], Mapping[str, object]]],
    expected: Set[Tuple[int, int, str, str, int, int]],
) -> Dict[str, object]:
    if not expected or set(records) != set(ARMS) or \
            any(set(records[arm]) != expected for arm in ARMS):
        die("evaluation cube does not exactly cover its sealed axes")
    axes: Dict[str, Dict[str, Dict[str, object]]] = {
        arm: {} for arm in ARMS}
    strata: Dict[str, Dict[str, Dict[str, object]]] = {
        arm: {} for arm in ARMS}
    arm_totals: Dict[str, Dict[str, object]] = {}
    weak_multiplicity: Dict[str, Dict[int, int]] = {arm: {} for arm in ARMS}
    weak_by_width: Dict[str, Dict[int, Set[int]]] = {
        arm: {width: set() for width in WIDTHS} for arm in ARMS}
    for arm in ARMS:
        total_trials = failures = errors = 0
        xors = Decimal(0)
        muladds = Decimal(0)
        for key, row in records[arm].items():
            K, root_index, schedule, loss, width, overhead = key
            row_trials = int(row["trials"])
            row_failures = int(row["rank_fail"])
            row_errors = int(row["error"])
            row_xors = Decimal(row["block_xors_mu"]) * row_trials
            row_muladds = Decimal(row["block_muladds_mu"]) * row_trials
            total_trials += row_trials
            failures += row_failures
            errors += row_errors
            xors += row_xors
            muladds += row_muladds
            if row_failures or row_errors:
                weak_multiplicity[arm][K] = weak_multiplicity[arm].get(K, 0) + 1
                weak_by_width[arm][width].add(K)
            for name, value in (
                    ("root", root_index), ("width", width),
                    ("schedule", schedule),
                    ("loss", loss), ("overhead", overhead)):
                axis = axes[arm].setdefault(_axis_key(name, value), {
                    "trials": 0, "rank_fail": 0, "error": 0,
                    "block_xors_sum": Decimal(0),
                    "block_muladds_sum": Decimal(0)})
                axis["trials"] = int(axis["trials"]) + row_trials
                axis["rank_fail"] = int(axis["rank_fail"]) + row_failures
                axis["error"] = int(axis["error"]) + row_errors
                axis["block_xors_sum"] = Decimal(
                    axis["block_xors_sum"]) + row_xors
                axis["block_muladds_sum"] = Decimal(
                    axis["block_muladds_sum"]) + row_muladds
            stratum_name = (
                "width={}|schedule={}|loss={}|overhead={}".format(
                    width, schedule, loss, overhead))
            stratum = strata[arm].setdefault(stratum_name, {
                "trials": 0, "rank_fail": 0, "error": 0,
                "block_xors_sum": Decimal(0),
                "block_muladds_sum": Decimal(0)})
            stratum["trials"] = int(stratum["trials"]) + row_trials
            stratum["rank_fail"] = int(stratum["rank_fail"]) + row_failures
            stratum["error"] = int(stratum["error"]) + row_errors
            stratum["block_xors_sum"] = Decimal(
                stratum["block_xors_sum"]) + row_xors
            stratum["block_muladds_sum"] = Decimal(
                stratum["block_muladds_sum"]) + row_muladds
        arm_totals[arm] = {
            "trials": total_trials, "rank_fail": failures, "error": errors,
            "failure_rate": str(Decimal(failures + errors) / total_trials),
            "weak_K": len(weak_multiplicity[arm]),
            "maximum_weak_K_multiplicity": max(
                weak_multiplicity[arm].values(), default=0),
            "block_xors_sum": str(xors), "block_muladds_sum": str(muladds),
        }
        for collection in (axes[arm], strata[arm]):
            for aggregate in collection.values():
                aggregate["block_xors_sum"] = str(
                    aggregate["block_xors_sum"])
                aggregate["block_muladds_sum"] = str(
                    aggregate["block_muladds_sum"])
    quadrants = {
        "both_success": 0, "repair": 0, "introduction": 0, "both_fail": 0,
    }
    repairs_by_overhead = {str(value): 0 for value in EVALUATION_OVERHEADS}
    introductions_by_overhead = {str(value): 0 for value in EVALUATION_OVERHEADS}
    repairs_by_root = {
        str(value): 0 for value in range(EVALUATION_ROOT_COUNT)}
    introductions_by_root = {
        str(value): 0 for value in range(EVALUATION_ROOT_COUNT)}
    paired_by_stratum: Dict[str, Dict[str, int]] = {}
    for key in sorted(expected):
        control = records["control"][key]
        candidate = records["candidate"][key]
        if control["trials"] != candidate["trials"]:
            die("paired evaluation trials differ")
        if int(control["error"]) or int(candidate["error"]):
            die("evaluation contains a codec error")
        c_fail = set(control["failures"])
        n_fail = set(candidate["failures"])
        both = len(c_fail.intersection(n_fail))
        repair = len(c_fail - n_fail)
        introduction = len(n_fail - c_fail)
        trials = int(control["trials"])
        quadrants["both_fail"] += both
        quadrants["repair"] += repair
        quadrants["introduction"] += introduction
        quadrants["both_success"] += trials - len(c_fail.union(n_fail))
        overhead = str(key[-1])
        repairs_by_overhead[overhead] += repair
        introductions_by_overhead[overhead] += introduction
        root_index = str(key[1])
        repairs_by_root[root_index] += repair
        introductions_by_root[root_index] += introduction
        stratum_name = "width={}|schedule={}|loss={}|overhead={}".format(
            key[-2], key[2], key[3], key[-1])
        paired = paired_by_stratum.setdefault(stratum_name, {
            "both_success": 0, "repair": 0,
            "introduction": 0, "both_fail": 0})
        paired["both_fail"] += both
        paired["repair"] += repair
        paired["introduction"] += introduction
        paired["both_success"] += trials - len(c_fail.union(n_fail))
    aggregate_nonworse = all(
        axes["candidate"][_axis_key("overhead", overhead)]["rank_fail"] <=
        axes["control"][_axis_key("overhead", overhead)]["rank_fail"]
        for overhead in EVALUATION_OVERHEADS)
    return {
        "arms": arm_totals, "axis_totals": axes,
        "width_schedule_loss_overhead_totals": strata,
        "weak_K_multiplicity": {
            arm: {str(K): weak_multiplicity[arm][K]
                  for K in sorted(weak_multiplicity[arm])}
            for arm in ARMS},
        "weak_K_by_width": {
            arm: {str(width): sorted(weak_by_width[arm][width])
                  for width in WIDTHS} for arm in ARMS},
        "paired_outcomes": quadrants,
        "paired_outcomes_by_width_schedule_loss_overhead": paired_by_stratum,
        "repairs_by_overhead": repairs_by_overhead,
        "introductions_by_overhead": introductions_by_overhead,
        "repairs_by_root": repairs_by_root,
        "introductions_by_root": introductions_by_root,
        "net_failure_change": (
            int(arm_totals["candidate"]["rank_fail"]) -
            int(arm_totals["control"]["rank_fail"])),
        "raw_recovery_gate": {
            "codec_errors_zero": all(
                int(arm_totals[arm]["error"]) == 0 for arm in ARMS),
            "aggregate_each_overhead_nonworse": aggregate_nonworse,
            "passed": aggregate_nonworse and all(
                int(arm_totals[arm]["error"]) == 0 for arm in ARMS),
            "per_width_schedule_loss_is_report_only": True,
            "weak_K_and_introductions_are_report_only_before_seed_fixes": True,
        },
    }


def reduce_evaluation(root: Path) -> Dict[str, object]:
    config = verify_config(root)
    cohort = verify_cohort_record(root, config, verify_discovery=True)
    jobs = load_phase_jobs(root, "evaluation", config)
    records, receipt_ledger = collect_phase_records(
        root, "evaluation", config, jobs)
    expected = {
        (K, root_index, schedule, loss, width, overhead)
        for K in cohort["cohort"]
        for root_index in range(EVALUATION_ROOT_COUNT)
        for schedule in SCHEDULES for loss in EVALUATION_LOSSES
        for width in WIDTHS for overhead in EVALUATION_OVERHEADS
    }
    summary = summarize_evaluation_records(records, expected)
    result = sealed_record(SCHEMA + ".result", {
        "config_sha256": sha256_file(root / "config.json"),
        "cohort_sha256": sha256_file(root / "cohort.json"),
        "evaluation_manifest_sha256": sha256_file(
            root / "evaluation" / "tasks.jsonl"),
        "evaluation_receipt_ledger_sha256": sha256_bytes(
            canonical_json(receipt_ledger)),
        **summary,
        "promotion_ready": False,
        "promotion_blocker": (
            "recovery-only evidence; independent full-payload speed gate and "
            "architecture decision are required"),
        "bb64_prior_evidence_scope": (
            "the degree-balanced all-K result used only BlockBytes=64; this "
            "result is the required natural-width cross-payload evidence"),
    })
    write_once(root / "evaluation" / "result.json", canonical_json(result))
    return result


def command_prepare(args: argparse.Namespace) -> None:
    root = resolve_direct_path(args.root, "campaign root")
    design_path = resolve_direct_path(args.design, "input design")
    controller = Path(__file__).resolve()
    if root.exists() or root.is_symlink():
        die("campaign root already exists")
    input_record, arms = load_input_design(design_path)
    staging = root.parent / (root.name + ".prepare.{}".format(os.getpid()))
    if staging.exists() or staging.is_symlink():
        die("campaign staging root already exists")
    try:
        staging.mkdir(parents=True)
        for directory in ("raw", "stderr", "receipts"):
            (staging / "discovery" / directory).mkdir(parents=True)
        config = prepare_config(staging, input_record, arms, controller)
        write_once(staging / "config.json", canonical_json(config))
        jobs = build_discovery_jobs(config)
        validate_jobs(jobs, "discovery", config)
        write_once(staging / "discovery" / "tasks.jsonl", manifest_bytes(jobs))
        audit = discovery_plan_record(staging, jobs)
        write_once(staging / "discovery" / "plan.json", canonical_json(audit))
        os.rename(str(staging), str(root))
    except Exception:
        shutil.rmtree(str(staging), ignore_errors=True)
        raise
    print(json.dumps({
        "status": "PREPARED_NOT_LAUNCHED", "root": str(root),
        "config_sha256": sha256_file(root / "config.json"),
        "discovery_manifest_sha256": sha256_file(
            root / "discovery" / "tasks.jsonl"),
        "jobs": len(jobs),
    }, indent=2, sort_keys=True))


def command_run(args: argparse.Namespace) -> None:
    root = resolve_direct_path(args.root, "campaign root")
    print(json.dumps(run_phase(
        root, args.phase, args.workers, args.execute), indent=2, sort_keys=True))


def command_seal_cohort(args: argparse.Namespace) -> None:
    cohort = seal_discovery_cohort(resolve_direct_path(args.root, "campaign root"))
    print(json.dumps({
        "status": "COHORT_SEALED_EVALUATION_NOT_LAUNCHED",
        "weak_union": len(cohort["weak_union"]),
        "controls": len(cohort["controls"]),
        "cohort": len(cohort["cohort"]),
        "cohort_sha256": sha256_bytes(canonical_json(cohort)),
    }, indent=2, sort_keys=True))


def command_reduce(args: argparse.Namespace) -> None:
    result = reduce_evaluation(resolve_direct_path(args.root, "campaign root"))
    print(json.dumps(result, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="freeze architecture and discovery; no launch")
    prepare.add_argument("--root", required=True)
    prepare.add_argument("--design", required=True)
    prepare.set_defaults(function=command_prepare)
    run = sub.add_parser("run", help="preflight or explicitly execute one frozen phase")
    run.add_argument("--root", required=True)
    run.add_argument("--phase", choices=("discovery", "evaluation"), required=True)
    run.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    run.add_argument("--execute", action="store_true")
    run.set_defaults(function=command_run)
    cohort = sub.add_parser(
        "seal-cohort", help="derive weak-K union and freeze evaluation ledger")
    cohort.add_argument("--root", required=True)
    cohort.set_defaults(function=command_seal_cohort)
    reduce_parser = sub.add_parser("reduce", help="strictly reduce completed evaluation")
    reduce_parser.add_argument("--root", required=True)
    reduce_parser.set_defaults(function=command_reduce)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        args.function(args)
    except CampaignError as exc:
        print("campaign error: {}".format(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
