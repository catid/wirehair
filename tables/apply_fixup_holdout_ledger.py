#!/usr/bin/env python3
"""Apply a validated exact-fixup ledger to new include files without clobbering.

The input include files must exactly match the ledger baseline.  Removal,
retention, and independently validated retunes are supported.  Any pending or
failed retune blocks the entire two-file publication.
"""

import argparse
import math
import re
import sys
from pathlib import Path

try:
    from . import atomic_publish
    from . import fixup_holdout_ledger as holdout
except ImportError:
    import atomic_publish
    import fixup_holdout_ledger as holdout


class ApplyError(ValueError):
    pass


_INTEGER_FIELDS = (
    "N", "PeelFixup", "DenseFixup", "BaselinePeelFixup",
    "BaselineDenseFixup", "BaseFail10", "BaseFail30", "FixedFail10",
    "FixedFail30", "BaseTrials10", "BaseTrials30", "FixedTrials10",
    "FixedTrials30", "MinimumTrials", "BaseSeed10", "BaseSeed30",
    "FixedSeed10", "FixedSeed30", "BlockBytes", "StartMode",
)
_FLOAT_FIELDS = (
    "BaseMean10", "BaseMean30", "FixedMean10", "FixedMean30",
    "MeanGate10", "MeanGate30",
)
_COMPLETE_DECISIONS = {
    "remove-exact-fixup",
    "retain-exact-fixup",
    "retune-exact-fixup",
}


def _passes(row, prefix):
    return (
        row[f"{prefix}Fail10"] == 0 and row[f"{prefix}Fail30"] == 0 and
        row[f"{prefix}Mean10"] <= row["MeanGate10"] and
        row[f"{prefix}Mean30"] <= row["MeanGate30"]
    )


def read_ledger(path):
    lines = Path(path).read_text(encoding="ascii").splitlines()
    if not lines or tuple(lines[0].split("\t")) != holdout.LEDGER_HEADER:
        raise ApplyError(f"{path}: unexpected ledger header")
    rows = []
    previous_n = 1
    methodology = None
    for line_number, line in enumerate(lines[1:], 2):
        if not line:
            raise ApplyError(f"{path}:{line_number}: blank ledger row")
        values = line.split("\t")
        if len(values) != len(holdout.LEDGER_HEADER):
            raise ApplyError(f"{path}:{line_number}: wrong ledger column count")
        row = dict(zip(holdout.LEDGER_HEADER, values))
        try:
            for field in _INTEGER_FIELDS:
                row[field] = int(row[field], 10)
            for field in _FLOAT_FIELDS:
                row[field] = float(row[field])
        except ValueError as exc:
            raise ApplyError(f"{path}:{line_number}: malformed numeric field") from exc
        n = row["N"]
        seeds = (
            row["PeelFixup"], row["DenseFixup"],
            row["BaselinePeelFixup"], row["BaselineDenseFixup"],
        )
        if (
            n <= previous_n or n > 64000 or
            any(seed < -1 or seed > 255 for seed in seeds) or
            any(not math.isfinite(row[field]) or row[field] < 0
                for field in _FLOAT_FIELDS) or
            any(row[field] < 0 for field in _INTEGER_FIELDS[5:]) or
            not row["Rationale"]
        ):
            raise ApplyError(f"{path}:{line_number}: invalid ledger row")
        if seeds == (-1, -1, -1, -1):
            raise ApplyError(f"{path}:{line_number}: ledger row has no exact fixup")
        previous_n = n
        if (
            row["MinimumTrials"] < 10000 or
            any(row[field] < row["MinimumTrials"] for field in (
                "BaseTrials10", "BaseTrials30", "FixedTrials10", "FixedTrials30"
            )) or
            row["MeanGate10"] > 0.05 or row["MeanGate30"] > 0.15
        ):
            raise ApplyError(
                f"{path}:{line_number}: holdout trials or mean gates are too weak"
            )
        if row["BlockBytes"] != 64 or row["StartMode"] != 0:
            raise ApplyError(
                f"{path}:{line_number}: holdouts must use bb=64 startMode=0"
            )
        if any(row[field] > (1 << 64) - 1 for field in (
            "BaseSeed10", "BaseSeed30", "FixedSeed10", "FixedSeed30"
        )):
            raise ApplyError(f"{path}:{line_number}: holdout seed exceeds uint64")
        if not re.fullmatch(r"[0-9]+(?:,[0-9]+)*", row["TrainingSeeds"]):
            raise ApplyError(f"{path}:{line_number}: malformed TrainingSeeds")
        training_seeds = {
            int(value, 10) for value in row["TrainingSeeds"].split(",")
        }
        if any(value > (1 << 64) - 1 for value in training_seeds):
            raise ApplyError(f"{path}:{line_number}: training seed exceeds uint64")
        try:
            holdout.validate_holdout_seeds(
                row["BaseSeed10"], row["BaseSeed30"],
                row["FixedSeed10"], row["FixedSeed30"], training_seeds,
            )
        except holdout.LedgerError as exc:
            raise ApplyError(f"{path}:{line_number}: {exc}") from exc
        row_methodology = (
            row["BaseTrials10"], row["BaseTrials30"],
            row["FixedTrials10"], row["FixedTrials30"],
            row["MinimumTrials"], row["MeanGate10"], row["MeanGate30"],
            row["BaseSeed10"], row["BaseSeed30"],
            row["FixedSeed10"], row["FixedSeed30"],
            row["BlockBytes"], row["StartMode"], row["TrainingSeeds"],
        )
        if methodology is None:
            methodology = row_methodology
        elif methodology != row_methodology:
            raise ApplyError(f"{path}:{line_number}: mixed holdout methodology")
        if (
            row["BaseFail10"] > row["BaseTrials10"] or
            row["BaseFail30"] > row["BaseTrials30"] or
            row["FixedFail10"] > row["FixedTrials10"] or
            row["FixedFail30"] > row["FixedTrials30"]
        ):
            raise ApplyError(f"{path}:{line_number}: failures exceed trials")

        candidate = row["PeelFixup"], row["DenseFixup"]
        baseline = row["BaselinePeelFixup"], row["BaselineDenseFixup"]
        base_good = _passes(row, "Base")
        fixed_good = _passes(row, "Fixed")
        decision = row["Decision"]
        consistent = (
            (decision == "remove-exact-fixup" and base_good) or
            (decision == "retain-exact-fixup" and not base_good and
             fixed_good and candidate == baseline) or
            (decision == "retune-exact-fixup" and not base_good and
             fixed_good and candidate != baseline and candidate != (-1, -1))
        )
        if decision not in _COMPLETE_DECISIONS:
            raise ApplyError(
                f"{path}:{line_number}: unresolved decision {decision!r} blocks application"
            )
        if not consistent:
            raise ApplyError(
                f"{path}:{line_number}: decision contradicts holdout evidence"
            )
        rows.append(row)
    if not rows:
        raise ApplyError(f"{path}: empty ledger")
    return rows


def apply_rows(rows, baseline_peel, baseline_dense):
    ledger_baseline_peel = {
        row["N"]: row["BaselinePeelFixup"]
        for row in rows if row["BaselinePeelFixup"] >= 0
    }
    ledger_baseline_dense = {
        row["N"]: row["BaselineDenseFixup"]
        for row in rows if row["BaselineDenseFixup"] >= 0
    }
    if baseline_peel != ledger_baseline_peel:
        raise ApplyError("peel input does not exactly match ledger baseline")
    if baseline_dense != ledger_baseline_dense:
        raise ApplyError("dense input does not exactly match ledger baseline")

    output_peel = dict(baseline_peel)
    output_dense = dict(baseline_dense)
    for row in rows:
        n = row["N"]
        if row["Decision"] == "remove-exact-fixup":
            output_peel.pop(n, None)
            output_dense.pop(n, None)
            continue
        for output, field in (
            (output_peel, "PeelFixup"), (output_dense, "DenseFixup")
        ):
            seed = row[field]
            if seed < 0:
                output.pop(n, None)
            else:
                output[n] = seed
    return output_peel, output_dense


def render_fixups(values, kind):
    lines = [
        f"// Exact-N {kind}-seed corrections validated by independent holdouts.\n",
        "// Sorted ascending by N (binary search). Format: { N, seed }\n",
    ]
    lines.extend(
        f"    {{ {n:5d}, {seed:3d} }},\n" for n, seed in sorted(values.items())
    )
    return "".join(lines)


def write_pair_no_replace(peel_path, peel_content, dense_path, dense_content):
    try:
        atomic_publish.publish_files_no_replace([
            (peel_path, peel_content),
            (dense_path, dense_content),
        ])
    except atomic_publish.AtomicPublishError as exc:
        raise ApplyError(str(exc)) from exc


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger", required=True)
    parser.add_argument("--peel-input", required=True)
    parser.add_argument("--dense-input", required=True)
    parser.add_argument("--peel-output", required=True)
    parser.add_argument("--dense-output", required=True)
    args = parser.parse_args(argv)
    try:
        rows = read_ledger(args.ledger)
        baseline_peel = holdout.read_fixups(args.peel_input)
        baseline_dense = holdout.read_fixups(args.dense_input)
        output_peel, output_dense = apply_rows(rows, baseline_peel, baseline_dense)
        write_pair_no_replace(
            args.peel_output,
            render_fixups(output_peel, "peel"),
            args.dense_output,
            render_fixups(output_dense, "dense"),
        )
    except (OSError, ApplyError, holdout.LedgerError) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(
        f"wrote peel={len(output_peel)} dense={len(output_dense)} exact fixups",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
