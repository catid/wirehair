#!/usr/bin/env python3
"""Build an auditable exact-N seed-fixup decision ledger.

A fixup is accepted only when both independent holdouts report zero codec
failures and mean overhead at or below 0.05 (10% loss) and 0.15 (30% loss).
Both the fixup-free and fixup-enabled measurements, their trial counts, and
the gates are retained in the durable TSV output.
"""

import argparse
import math
import re
import sys
from pathlib import Path

try:
    from . import atomic_publish
except ImportError:
    import atomic_publish


class LedgerError(ValueError):
    pass


MASK64 = (1 << 64) - 1


LEDGER_HEADER = (
    "N", "PeelFixup", "DenseFixup", "BaselinePeelFixup",
    "BaselineDenseFixup", "BaseMean10", "BaseMean30", "BaseFail10",
    "BaseFail30", "FixedMean10", "FixedMean30", "FixedFail10",
    "FixedFail30", "BaseTrials10", "BaseTrials30", "FixedTrials10",
    "FixedTrials30", "MinimumTrials", "MeanGate10", "MeanGate30",
    "BaseSeed10", "BaseSeed30", "FixedSeed10", "FixedSeed30",
    "BlockBytes", "StartMode", "TrainingSeeds", "Decision", "Rationale",
)


def read_fixups(path):
    source = Path(path).read_text(encoding="ascii")
    source = re.sub(r"//[^\n]*", "", source)
    pairs = [
        (int(n, 10), int(seed, 10))
        for n, seed in re.findall(r"\{\s*([0-9]+)\s*,\s*([0-9]+)\s*\}\s*,", source)
    ]
    cleaned = re.sub(r"\{\s*[0-9]+\s*,\s*[0-9]+\s*\}\s*,|\s", "", source)
    if cleaned:
        raise LedgerError(f"{path}: unsupported fixup syntax: {cleaned[:40]!r}")
    previous = 1
    values = {}
    for n, seed in pairs:
        if n <= previous or n > 64000 or seed > 255:
            raise LedgerError(f"{path}: unsorted or out-of-range fixup N={n} seed={seed}")
        previous = n
        values[n] = seed
    return values


def read_ohead(
    path, expected_loss, minimum_trials, expected_bb=64, expected_start_mode=0
):
    lines = Path(path).read_text(encoding="ascii").splitlines()
    metadata = None
    header_seen = False
    rows = {}
    for line_number, line in enumerate(lines, 1):
        if line.startswith("# ohead:"):
            if metadata is not None:
                raise LedgerError(f"{path}:{line_number}: duplicate ohead metadata")
            trials = re.findall(r"trials/N=([0-9]+)", line)
            loss = re.findall(r"loss=([0-9]+(?:\.[0-9]+)?)", line)
            fixups = re.findall(r"seed_fixups=(enabled|disabled)", line)
            count = re.findall(r"Ns=([0-9]+)", line)
            block_bytes = re.findall(r"bb=([0-9]+)", line)
            start_mode = re.findall(r"startMode=([0-9]+)", line)
            seed = re.findall(r"seed=(0x[0-9a-fA-F]+|[0-9]+)", line)
            fields = (
                trials, loss, fixups, count, block_bytes, start_mode, seed,
            )
            if any(len(values) != 1 for values in fields):
                raise LedgerError(
                    f"{path}:{line_number}: malformed ohead metadata; "
                    "exactly one trials, loss, seed, seed_fixups, Ns, bb, and "
                    "startMode field is required"
                )
            metadata = (
                int(trials[0]), float(loss[0]), fixups[0], int(count[0]),
                int(block_bytes[0]), int(start_mode[0]), int(seed[0], 0),
            )
            continue
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if fields == ["N", "mean", "p50", "p99", "p999", "max", "fail"]:
            if header_seen:
                raise LedgerError(f"{path}:{line_number}: duplicate ohead header")
            header_seen = True
            continue
        if not header_seen:
            raise LedgerError(f"{path}:{line_number}: ohead row precedes header")
        if len(fields) != 7:
            raise LedgerError(f"{path}:{line_number}: expected seven ohead columns")
        try:
            n = int(fields[0], 10)
            mean = float(fields[1])
            quantiles = [int(value, 10) for value in fields[2:6]]
            fail = int(fields[6], 10)
        except ValueError as exc:
            raise LedgerError(f"{path}:{line_number}: malformed ohead row") from exc
        if (
            n in rows or n < 2 or n > 64000 or
            not math.isfinite(mean) or mean < 0 or fail < 0 or
            any(value < 0 for value in quantiles) or
            quantiles != sorted(quantiles) or mean > quantiles[-1]
        ):
            raise LedgerError(f"{path}:{line_number}: invalid or duplicate ohead row")
        rows[n] = (mean, fail)
    if metadata is None or not header_seen or not rows:
        raise LedgerError(f"{path}: missing ohead metadata, header, or rows")
    trials, loss, policy, expected_count, block_bytes, start_mode, seed = metadata
    if seed > MASK64:
        raise LedgerError(f"{path}: holdout seed exceeds uint64")
    if trials < minimum_trials or abs(loss - expected_loss) > 1e-12:
        raise LedgerError(
            f"{path}: expected loss={expected_loss:.2f}, trials>={minimum_trials}; "
            f"got loss={loss:.2f}, trials={trials}"
        )
    if expected_count != len(rows):
        raise LedgerError(
            f"{path}: metadata Ns={expected_count} does not match {len(rows)} rows"
        )
    if any(failures > trials for _, failures in rows.values()):
        raise LedgerError(f"{path}: failure count exceeds trials/N={trials}")
    if block_bytes != expected_bb or block_bytes <= 0:
        raise LedgerError(
            f"{path}: expected real-byte bb={expected_bb}; got bb={block_bytes}"
        )
    if start_mode != expected_start_mode or start_mode != 0:
        raise LedgerError(
            f"{path}: loss holdout must use startMode=0; got {start_mode}"
        )
    return rows, policy, trials, seed, block_bytes, start_mode


def read_training_seeds(paths):
    seeds = set()
    for path in paths:
        found = False
        for line_number, line in enumerate(
            Path(path).read_text(encoding="ascii").splitlines(), 1
        ):
            match = re.fullmatch(r"(seed|seeds)=([0-9]+(?:,[0-9]+)*)", line)
            if not match:
                continue
            found = True
            for value in match.group(2).split(","):
                seed = int(value, 10)
                if seed > MASK64:
                    raise LedgerError(
                        f"{path}:{line_number}: training seed exceeds uint64"
                    )
                seeds.add(seed)
        if not found:
            raise LedgerError(f"{path}: no seed=/seeds= campaign parameter")
    if not seeds:
        raise LedgerError("no training seeds")
    return seeds


def validate_holdout_seeds(
    base_seed10, base_seed30, fixed_seed10, fixed_seed30, training_seeds
):
    if base_seed10 != fixed_seed10 or base_seed30 != fixed_seed30:
        raise LedgerError("base/fixed holdouts must pair the same seed within each loss")
    if base_seed10 == base_seed30:
        raise LedgerError("loss10 and loss30 holdouts must use distinct seeds")
    reused = {base_seed10, base_seed30} & set(training_seeds)
    if reused:
        raise LedgerError(
            f"holdout seeds overlap generator/training seeds: {sorted(reused)}"
        )


def _fixup_pair(peel, dense, n):
    return peel.get(n, -1), dense.get(n, -1)


def _passes(rows10, rows30, n, good10, good30):
    mean10, failures10 = rows10[n]
    mean30, failures30 = rows30[n]
    return (
        failures10 == 0 and failures30 == 0 and
        mean10 <= good10 and mean30 <= good30
    )


def build_ledger(
    peel_fixups,
    dense_fixups,
    base10,
    base30,
    fixed10=None,
    fixed30=None,
    good10=0.05,
    good30=0.15,
    *,
    baseline_peel_fixups=None,
    baseline_dense_fixups=None,
    base_trials10=0,
    base_trials30=0,
    fixed_trials10=0,
    fixed_trials30=0,
    minimum_trials=0,
    base_seed10=0,
    base_seed30=0,
    fixed_seed10=0,
    fixed_seed30=0,
    block_bytes=64,
    start_mode=0,
    training_seeds=(),
):
    if baseline_peel_fixups is None:
        baseline_peel_fixups = peel_fixups
    if baseline_dense_fixups is None:
        baseline_dense_fixups = dense_fixups
    expected = (
        set(peel_fixups) | set(dense_fixups) |
        set(baseline_peel_fixups) | set(baseline_dense_fixups)
    )
    if not expected:
        raise LedgerError("exact fixup union is empty")
    for label, rows in (("base10", base10), ("base30", base30)):
        if set(rows) != expected:
            missing = sorted(expected - set(rows))
            extra = sorted(set(rows) - expected)
            raise LedgerError(
                f"{label} coverage mismatch: missing={missing[:8]} extra={extra[:8]}"
            )
    if (fixed10 is None) != (fixed30 is None):
        raise LedgerError("enabled-fixup holdouts must be provided as a pair")
    if fixed10 is not None:
        for label, rows in (("fixed10", fixed10), ("fixed30", fixed30)):
            if set(rows) != expected:
                raise LedgerError(f"{label} coverage does not match exact fixup union")

    output = []
    for n in sorted(expected):
        candidate = _fixup_pair(peel_fixups, dense_fixups, n)
        baseline = _fixup_pair(baseline_peel_fixups, baseline_dense_fixups, n)
        base_good = _passes(base10, base30, n, good10, good30)
        fixed_good = fixed10 is not None and _passes(
            fixed10, fixed30, n, good10, good30
        )
        if base_good:
            decision = "remove-exact-fixup"
            rationale = "fixup-free holdouts pass both zero-failure mean gates"
        elif fixed10 is None:
            decision = "retain-pending-enabled-holdout"
            rationale = "fixup-free holdout misses a gate; enabled holdouts are absent"
        elif not fixed_good:
            decision = "retune-required"
            rationale = "fixup-free and fixup-enabled holdouts both miss a gate"
        elif candidate == baseline:
            decision = "retain-exact-fixup"
            rationale = "base misses a gate; unchanged exact fixup passes both gates"
        elif candidate != (-1, -1):
            decision = "retune-exact-fixup"
            rationale = "base misses a gate; changed exact fixup passes both gates"
        else:
            decision = "retune-required"
            rationale = "base misses a gate and the candidate removed all exact fixups"

        base_mean10, base_fail10 = base10[n]
        base_mean30, base_fail30 = base30[n]
        if fixed10 is None:
            fixed_mean10 = fixed_mean30 = None
            fixed_fail10 = fixed_fail30 = None
        else:
            fixed_mean10, fixed_fail10 = fixed10[n]
            fixed_mean30, fixed_fail30 = fixed30[n]
        output.append({
            "N": n,
            "PeelFixup": candidate[0],
            "DenseFixup": candidate[1],
            "BaselinePeelFixup": baseline[0],
            "BaselineDenseFixup": baseline[1],
            "BaseMean10": base_mean10,
            "BaseMean30": base_mean30,
            "BaseFail10": base_fail10,
            "BaseFail30": base_fail30,
            "FixedMean10": fixed_mean10,
            "FixedMean30": fixed_mean30,
            "FixedFail10": fixed_fail10,
            "FixedFail30": fixed_fail30,
            "BaseTrials10": base_trials10,
            "BaseTrials30": base_trials30,
            "FixedTrials10": fixed_trials10,
            "FixedTrials30": fixed_trials30,
            "MinimumTrials": minimum_trials,
            "MeanGate10": good10,
            "MeanGate30": good30,
            "BaseSeed10": base_seed10,
            "BaseSeed30": base_seed30,
            "FixedSeed10": fixed_seed10,
            "FixedSeed30": fixed_seed30,
            "BlockBytes": block_bytes,
            "StartMode": start_mode,
            "TrainingSeeds": ",".join(str(value) for value in sorted(training_seeds)),
            "Decision": decision,
            "Rationale": rationale,
        })
    return output


def _format_value(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return format(value, ".17g")
    return str(value)


def write_atomic(path, rows):
    output = ["\t".join(LEDGER_HEADER) + "\n"]
    for row in rows:
        values = [_format_value(row[column]) for column in LEDGER_HEADER]
        if any("\t" in value or "\n" in value for value in values):
            raise LedgerError("ledger value contains a tab or newline")
        output.append("\t".join(values) + "\n")
    try:
        atomic_publish.publish_files_no_replace([(path, "".join(output))])
    except atomic_publish.AtomicPublishError as exc:
        raise LedgerError(str(exc)) from exc


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Require zero failures and mean overhead <=0.05/0.15 at loss "
            "0.10/0.30 for both explicit fixup policies."
        )
    )
    parser.add_argument("--peel-fixups", required=True)
    parser.add_argument("--dense-fixups", required=True)
    parser.add_argument("--baseline-peel-fixups")
    parser.add_argument("--baseline-dense-fixups")
    parser.add_argument("--base-loss10", required=True)
    parser.add_argument("--base-loss30", required=True)
    parser.add_argument("--fixed-loss10", required=True)
    parser.add_argument("--fixed-loss30", required=True)
    parser.add_argument("--minimum-trials", type=int, default=10000)
    parser.add_argument("--block-bytes", type=int, default=64)
    parser.add_argument(
        "--training-manifest", action="append", required=True,
        help="campaign manifest containing canonical seed=/seeds= parameters",
    )
    parser.add_argument("--good10", type=float, default=0.05)
    parser.add_argument("--good30", type=float, default=0.15)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        if (
            args.minimum_trials < 1 or
            args.block_bytes < 1 or
            not math.isfinite(args.good10) or not math.isfinite(args.good30) or
            not 0 <= args.good10 <= args.good30
        ):
            raise LedgerError("invalid trial/threshold parameters")
        if bool(args.baseline_peel_fixups) != bool(args.baseline_dense_fixups):
            raise LedgerError("baseline fixup files must be supplied as a pair")
        peel = read_fixups(args.peel_fixups)
        dense = read_fixups(args.dense_fixups)
        baseline_peel = (
            read_fixups(args.baseline_peel_fixups)
            if args.baseline_peel_fixups else peel
        )
        baseline_dense = (
            read_fixups(args.baseline_dense_fixups)
            if args.baseline_dense_fixups else dense
        )
        training_seeds = read_training_seeds(args.training_manifest)
        base10, policy10, base_trials10, base_seed10, bb10, mode10 = read_ohead(
            args.base_loss10, 0.10, args.minimum_trials, args.block_bytes, 0
        )
        base30, policy30, base_trials30, base_seed30, bb30, mode30 = read_ohead(
            args.base_loss30, 0.30, args.minimum_trials, args.block_bytes, 0
        )
        if policy10 != "disabled" or policy30 != "disabled":
            raise LedgerError("base holdouts must explicitly use seed fixups disabled")
        fixed10, fixed_policy10, fixed_trials10, fixed_seed10, fbb10, fmode10 = read_ohead(
            args.fixed_loss10, 0.10, args.minimum_trials, args.block_bytes, 0
        )
        fixed30, fixed_policy30, fixed_trials30, fixed_seed30, fbb30, fmode30 = read_ohead(
            args.fixed_loss30, 0.30, args.minimum_trials, args.block_bytes, 0
        )
        if fixed_policy10 != "enabled" or fixed_policy30 != "enabled":
            raise LedgerError("fixed holdouts must explicitly use seed fixups enabled")
        validate_holdout_seeds(
            base_seed10, base_seed30, fixed_seed10, fixed_seed30,
            training_seeds,
        )
        if {bb10, bb30, fbb10, fbb30} != {args.block_bytes} or \
                {mode10, mode30, fmode10, fmode30} != {0}:
            raise LedgerError("holdout bb/startMode methodology mismatch")
        ledger = build_ledger(
            peel,
            dense,
            base10,
            base30,
            fixed10,
            fixed30,
            args.good10,
            args.good30,
            baseline_peel_fixups=baseline_peel,
            baseline_dense_fixups=baseline_dense,
            base_trials10=base_trials10,
            base_trials30=base_trials30,
            fixed_trials10=fixed_trials10,
            fixed_trials30=fixed_trials30,
            minimum_trials=args.minimum_trials,
            base_seed10=base_seed10,
            base_seed30=base_seed30,
            fixed_seed10=fixed_seed10,
            fixed_seed30=fixed_seed30,
            block_bytes=args.block_bytes,
            start_mode=0,
            training_seeds=training_seeds,
        )
        write_atomic(args.output, ledger)
    except (OSError, LedgerError) as exc:
        print(exc, file=sys.stderr)
        return 1
    counts = {}
    for row in ledger:
        counts[row["Decision"]] = counts.get(row["Decision"], 0) + 1
    print(f"wrote {len(ledger)} fixup decisions: {counts}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
