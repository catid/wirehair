#!/usr/bin/env python3
"""rank_total.py -- heavy-aware total block-op ranking for precode_sim CSVs.

For every precode_sim result row and every block size bb, computes

  total_block_ops(bb) = recv_xors_per_packet * (K + oh)
                      + precode_gen_xors_mu
                      + sparse_solve_xors_mu
                      + backsub_xors_mu
                      + (ge_real_rowops_mu if present else ge_block_xors_mu)
                      + ratio(bb) * heavy_muladds_mu

where ratio(bb) = muladd_us(bb) / xor_us(bb) from the two calibration CSVs
(experiments/peeling/results/{muladd,cold_xor}_calibration.csv).  All terms
except the heavy one are block-XOR counts; the heavy term converts block-unit
GF(256) muladds into XOR-equivalents at the measured per-block-size ratio.
--pessimistic doubles the heavy term.

Failure is re-scored for alternative heavy counts H' in {6, 8, 12, 16} from
the def_pdf column (def > H' = failure).  def_pdf is H-independent by
construction (it is recorded before the heavy MDS patch is applied); the
script asserts this per row by checking that re-scoring at the row's modeled
H reproduces the recorded fail_rate.  Cost columns always reflect the row's
modeled H; rows where a fail@H' column differs from the modeled H are not
cost-adjusted (the modeled-H column is marked with '*').  Rows with nonzero
runaway_rate are sorted after finite-cost rows because their cost means are
conditional on completed trials only.

Output: a ranked table (ascending total_block_ops) per (K, oh, bb), plus the
muladd:xor crossover ratio at which the ldpcdense H=12 scheme stops beating
the dense H=6 baseline (totals are linear in the ratio, so the crossover is
block-size independent).

Examples:
  python3 rank_total.py results/hybrid_K1000.csv \
      --xor-csv ../peeling/results/cold_xor_calibration.csv \
      --muladd-csv ../peeling/results/muladd_calibration.csv
  python3 rank_total.py --self-test
"""

import argparse
import csv
import os
import sys
import tempfile

H_PRIMES = (6, 8, 12, 16)
DEFAULT_BB = "1280,102400,1048576"
PDF_TOLERANCE = 1e-4


def parse_bb_list(text):
    """Parse a non-empty CSV of positive decimal block sizes."""
    if text is None or text == "":
        raise ValueError("--bb must be a non-empty CSV of positive integers")
    out = []
    for part in text.split(","):
        if not part or not part.isdigit():
            raise ValueError("--bb contains an invalid token: %r" % part)
        value = int(part)
        if value <= 0:
            raise ValueError("--bb values must be positive")
        out.append(value)
    return out


def parse_def_pdf(text):
    """'6:0.983250|7:0.016000' -> {6: 0.98325, 7: 0.016}"""
    pdf = {}
    for part in text.split("|"):
        if not part:
            continue
        d, p = part.split(":")
        pdf[int(d)] = float(p)
    return pdf


def fail_at(pdf, h):
    """P(def > h) from a def pdf (clamped: 6-dp pdf entries can sum to >1)."""
    return min(1.0, sum(p for d, p in pdf.items() if d > h))


def load_us_per_op(path, label):
    """Calibration CSV -> {block_bytes: us_per_op}, or None if unusable.

    Accepts any header with a block_bytes column and a us-per-op column.
    Column selection is label-aware because muladd_calibration.csv carries
    BOTH xor and muladd timings (xor_median_us_per_op,
    muladd_median_us_per_op); cold_xor_calibration.csv has median_us_per_xor.
    """
    if not path:
        return None
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        sys.stderr.write(
            "rank_total: %s calibration %s missing/empty; skipping it\n"
            % (label, path))
        return None
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "block_bytes" not in reader.fieldnames:
            sys.stderr.write(
                "rank_total: %s calibration %s has no block_bytes column\n"
                % (label, path))
            return None
        us_col = None
        for name in reader.fieldnames:
            if "us_per" in name and label in name:
                us_col = name
                break
        if us_col is None:
            for name in reader.fieldnames:
                if name.startswith("median_us_per"):
                    us_col = name
                    break
        if us_col is None:
            for name in reader.fieldnames:
                if "us_per" in name:
                    us_col = name
                    break
        if us_col is None:
            sys.stderr.write(
                "rank_total: %s calibration %s has no us-per-op column\n"
                % (label, path))
            return None
        table = {}
        for row in reader:
            try:
                table[int(row["block_bytes"])] = float(row[us_col])
            except (KeyError, ValueError):
                continue
        return table or None


def build_ratios(bb_list, xor_us, muladd_us, assume_ratio):
    """{bb: muladd:xor ratio}; falls back to assume_ratio with a warning."""
    ratios = {}
    for bb in bb_list:
        if (xor_us and muladd_us and bb in xor_us and bb in muladd_us
                and xor_us[bb] > 0.0):
            ratios[bb] = muladd_us[bb] / xor_us[bb]
        else:
            sys.stderr.write(
                "rank_total: no calibration for bb=%d; assuming "
                "muladd:xor ratio %.2f\n" % (bb, assume_ratio))
            ratios[bb] = assume_ratio
    return ratios


def load_rows(paths):
    """Parse precode_sim CSVs into scored row dicts (base/heavy split)."""
    rows = []
    for path in paths:
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            fields = set(reader.fieldnames or ())
            required = {"K", "scheme", "H", "oh", "fail_rate", "def_pdf",
                        "inact_mu", "precode_gen_xors_mu",
                        "sparse_solve_xors_mu", "backsub_xors_mu",
                        "ge_block_xors_mu"}
            missing = required - fields
            if missing:
                raise ValueError("%s: missing columns: %s"
                                 % (path, ",".join(sorted(missing))))
            for raw in reader:
                rows.append(score_row(raw, fields, path))
    return rows


def score_row(raw, fields, path):
    k = int(raw["K"])
    oh = int(raw["oh"])
    h = int(raw["H"])
    inact_mu = float(raw["inact_mu"])

    # Receive-side term: per-packet column scaled by packet count, with a
    # derivation fallback for older CSV layouts.
    if "recv_xors_per_packet" in fields and raw.get("recv_xors_per_packet"):
        recv = float(raw["recv_xors_per_packet"]) * (k + oh)
    elif "recv_xors_mu" in fields and raw.get("recv_xors_mu"):
        recv = float(raw["recv_xors_mu"])
    else:
        sys.stderr.write(
            "rank_total: %s: no receive-XOR column; recv term = 0\n" % path)
        recv = 0.0

    # GE term: real replay row ops where present, else the R^2/2 estimate
    if "ge_real_rowops_mu" in fields and raw.get("ge_real_rowops_mu"):
        ge = float(raw["ge_real_rowops_mu"])
        ge_src = "real"
    else:
        ge = float(raw["ge_block_xors_mu"])
        ge_src = "est"

    # Heavy term: recorded column where present, else derive H*inact + H^2
    if "heavy_muladds_mu" in fields and raw.get("heavy_muladds_mu"):
        heavy = float(raw["heavy_muladds_mu"])
    else:
        heavy = h * inact_mu + h * h

    pdf = parse_def_pdf(raw["def_pdf"])
    mass = sum(pdf.values())
    if abs(mass - 1.0) > 1e-3:
        raise ValueError("%s: def_pdf mass %.6f != 1 for %s K=%d oh=%d"
                         % (path, mass, raw["scheme"], k, oh))
    fail_rate = float(raw["fail_rate"])
    runaway_rate = float(raw.get("runaway_rate") or 0.0)
    rescored = fail_at(pdf, h)
    # def_pdf must be H-independent by construction: re-scoring at the row's
    # own modeled H has to reproduce the recorded fail_rate.
    if abs(rescored - fail_rate) > PDF_TOLERANCE:
        raise AssertionError(
            "%s: def_pdf is not H-independent: fail@H=%d from pdf is %.6f "
            "but fail_rate is %.6f (%s K=%d oh=%d)"
            % (path, h, rescored, fail_rate, raw["scheme"], k, oh))

    base = (recv + float(raw["precode_gen_xors_mu"])
            + float(raw["sparse_solve_xors_mu"])
            + float(raw["backsub_xors_mu"]) + ge)
    return {
        "K": k, "oh": oh, "H": h, "scheme": raw["scheme"], "file": path,
        "base": base, "heavy": heavy, "ge_src": ge_src,
        "fail_rate": fail_rate, "runaway_rate": runaway_rate,
        "fails": {hp: fail_at(pdf, hp) for hp in H_PRIMES},
    }


def total_ops(row, ratio, pessimistic):
    if row.get("runaway_rate", 0.0) > 0.0:
        return float("inf")
    factor = 2.0 if pessimistic else 1.0
    return row["base"] + ratio * factor * row["heavy"]


def crossover_ratio(dense_row, ldpc_row, pessimistic):
    """muladd:xor ratio where ldpc total equals dense total.

    Returns (ratio, verdict).  The verdict states the direction: when ldpc
    has MORE heavy work it beats dense below the ratio; when ldpc has LESS
    heavy work but a heavier base it beats dense above the ratio.  ratio is
    None for the degenerate cases (never / always beats).
    """
    factor = 2.0 if pessimistic else 1.0
    dheavy = factor * (ldpc_row["heavy"] - dense_row["heavy"])
    dbase = dense_row["base"] - ldpc_row["base"]
    if dheavy > 0.0:
        r = dbase / dheavy
        if r <= 0.0:
            return None, "never beats (base cost already higher)"
        return r, "beats below this ratio"
    if dbase > 0.0:
        return None, "beats at every ratio (cheaper base and heavy)"
    if dheavy < 0.0 and dbase <= 0.0:
        # ldpc is heavier in base ops but lighter in heavy ops: the heavy
        # savings grow with the ratio, so ldpc wins ABOVE the crossover.
        return dbase / dheavy, "beats above this ratio"
    return None, "never beats (heavier base, no heavy advantage)"


def print_report(rows, bb_list, ratios, pessimistic, out=sys.stdout):
    mode = "pessimistic(2x heavy)" if pessimistic else "normal"
    files = sorted(set(r["file"] for r in rows))
    groups = {}
    for r in rows:
        groups.setdefault((r["K"], r["oh"]), []).append(r)

    for (k, oh) in sorted(groups):
        grp = groups[(k, oh)]
        for bb in bb_list:
            ratio = ratios[bb]
            ranked = sorted(grp, key=lambda r: (total_ops(r, ratio,
                                                          pessimistic),
                                                r["scheme"], r["file"]))
            out.write("== K=%d oh=%d bb=%d (muladd:xor=%.3f, %s) ==\n"
                      % (k, oh, bb, ratio, mode))
            hdr = " %2s  %-28s %3s " % ("#", "scheme", "H")
            hdr += "".join("%11s" % ("fail@H%d" % hp) for hp in H_PRIMES)
            hdr += "%14s %5s %12s %9s\n" % (
                "total_ops", "ge", "heavy_ops", "runaway")
            out.write(hdr)
            for i, r in enumerate(ranked):
                name = r["scheme"]
                if len(files) > 1:
                    name += " [%s]" % os.path.basename(r["file"])
                cells = ""
                for hp in H_PRIMES:
                    mark = "*" if hp == r["H"] else " "
                    cells += "  %.6f%s" % (r["fails"][hp], mark)
                factor = 2.0 if pessimistic else 1.0
                total = total_ops(r, ratio, pessimistic)
                total_text = "runaway" if total == float("inf") else "%.1f" % total
                out.write(" %2d  %-28s %3d %s%14s %5s %12.1f %9.6f\n"
                          % (i + 1, name, r["H"], cells,
                             total_text, r["ge_src"],
                             ratio * factor * r["heavy"],
                             r["runaway_rate"]))
            out.write("\n")

        # Crossover: block-size independent (totals are linear in the ratio)
        dense6 = [r for r in grp
                  if r["runaway_rate"] == 0.0 and r["H"] == 6
                  and (r["scheme"] == "dense"
                       or r["scheme"].startswith("dense_d")
                       or r["scheme"] == "dense_h6")]
        ldpc12 = [r for r in grp
                  if r["runaway_rate"] == 0.0 and r["H"] == 12
                  and r["scheme"].startswith("ldpcdense")]
        if dense6 and ldpc12:
            for lr in ldpc12:
                for dr in dense6:
                    r, verdict = crossover_ratio(dr, lr, pessimistic)
                    if r is not None:
                        out.write("crossover K=%d oh=%d: %s vs %s: "
                                  "muladd:xor=%.3f (%s)\n"
                                  % (k, oh, lr["scheme"], dr["scheme"],
                                     r, verdict))
                    else:
                        out.write("crossover K=%d oh=%d: %s vs %s: %s\n"
                                  % (k, oh, lr["scheme"], dr["scheme"],
                                     verdict))
        else:
            out.write("crossover K=%d oh=%d: skipped (need a dense H=6 and "
                      "an ldpcdense H=12 row)\n" % (k, oh))
        out.write("\n")

    out.write("notes: fail@H' re-scored from def_pdf (def > H' = failure); "
              "'*' marks each row's modeled H.  Cost columns always use the "
              "modeled H; fail@H' columns with H' != modeled H are NOT "
              "cost-adjusted.  Nonzero runaway_rate means --max-inact "
              "aborted some trials, so total_ops is reported as runaway and "
              "sorted after finite-cost rows.\n")


#------------------------------------------------------------------------------
# Self-test

SELFTEST_SIM = """K,scheme,D,H,oh,trials,fail_rate,fail_rate_noheavy,def_mu,def_max,def_pdf,inact_mu,inact_sd,inact_max,rank_mu,peeled_mu,residual_rows_mu,recv_xors_per_packet,precode_gen_xors_mu,sparse_solve_xors_mu,backsub_xors_mu,ge_block_xors_mu,ge_bitops_mu,heavy_muladds_mu,heavy_divs_mu,ge_real_bitops_mu,ge_real_rowops_mu,fill_in_mu,def_outside_w18_rate,def_band_w95,def_band_w99
100,dense,10,6,0,1000,0.100000,1.000000,1.0,13,0:0.900000|7:0.050000|13:0.050000,20.00,1.00,25,19.0,80.0,19.0,5.0000,100.0,200.0,300.0,50.0,999.0,156.0,6.0,400.0,40.0,2.00,0.000000,3,5
100,ldpcdense_s8_d4_h12,8,12,0,1000,0.050000,1.000000,1.0,20,0:0.950000|13:0.040000|20:0.010000,15.00,1.00,20,14.0,85.0,14.0,4.0000,80.0,150.0,250.0,45.0,888.0,324.0,12.0,300.0,30.0,1.00,0.000000,2,4
"""

# Same rows but without heavy/ge_real columns: tests derivation fallbacks
SELFTEST_SIM_OLD = """K,scheme,D,H,oh,trials,fail_rate,fail_rate_noheavy,def_mu,def_max,def_pdf,inact_mu,inact_sd,inact_max,rank_mu,peeled_mu,residual_rows_mu,recv_xors_per_packet,precode_gen_xors_mu,sparse_solve_xors_mu,backsub_xors_mu,ge_block_xors_mu,ge_bitops_mu
100,dense,10,6,0,1000,0.100000,1.000000,1.0,13,0:0.900000|7:0.050000|13:0.050000,20.00,1.00,25,19.0,80.0,19.0,5.0000,100.0,200.0,300.0,50.0,999.0
"""

SELFTEST_SIM_RUNAWAY = """K,scheme,D,H,oh,trials,fail_rate,fail_rate_noheavy,def_mu,def_max,def_pdf,inact_mu,inact_sd,inact_max,rank_mu,peeled_mu,residual_rows_mu,recv_xors_per_packet,precode_gen_xors_mu,sparse_solve_xors_mu,backsub_xors_mu,ge_block_xors_mu,ge_bitops_mu,heavy_muladds_mu,heavy_divs_mu,runaway_rate
100,ldpcdense_s8_d4,8,6,0,1000,1.000000,1.000000,0.0,0,999999:1.000000,0.00,0.00,0,0.0,0.0,0.0,0.0000,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.000000
"""

SELFTEST_XOR = """block_bytes,working_blocks,ops_per_repeat,repeats,median_total_us,median_us_per_xor,median_gib_per_s,checksum
1280,1,1,1,1.0,1.0,1.0,0x0
102400,1,1,1,1.0,2.0,1.0,0x0
"""

# Real muladd_calibration.csv layout: BOTH xor and muladd timing columns.
# The xor_median_us_per_op values are decoys: the parser must pick the
# muladd_median_us_per_op column for the "muladd" label.
SELFTEST_MULADD = """block_bytes,pool_blocks,ops_per_repeat,repeats,xor_median_us_per_op,xor_median_gib_per_s,muladd_median_us_per_op,muladd_median_gib_per_s,muladd_xor_ratio,checksum
1280,1,1,1,1.5,1.0,3.0,1.0,2.0,0x0
102400,1,1,1,1.5,1.0,8.0,1.0,5.33,0x0
"""


def approx(a, b, eps=1e-9):
    return abs(a - b) <= eps


def self_test():
    tmpdir = tempfile.mkdtemp(prefix="rank_total_selftest_")
    sim_path = os.path.join(tmpdir, "sim.csv")
    old_path = os.path.join(tmpdir, "sim_old.csv")
    xor_path = os.path.join(tmpdir, "xor.csv")
    mul_path = os.path.join(tmpdir, "muladd.csv")
    for path, text in ((sim_path, SELFTEST_SIM), (old_path, SELFTEST_SIM_OLD),
                       (xor_path, SELFTEST_XOR), (mul_path, SELFTEST_MULADD)):
        with open(path, "w") as f:
            f.write(text)
    runaway_path = os.path.join(tmpdir, "runaway.csv")
    with open(runaway_path, "w") as f:
        f.write(SELFTEST_SIM_RUNAWAY)

    xor_us = load_us_per_op(xor_path, "xor")
    mul_us = load_us_per_op(mul_path, "muladd")
    assert xor_us == {1280: 1.0, 102400: 2.0}, xor_us
    assert mul_us == {1280: 3.0, 102400: 8.0}, mul_us
    ratios = build_ratios([1280, 102400, 999], xor_us, mul_us, 5.0)
    assert approx(ratios[1280], 3.0) and approx(ratios[102400], 4.0)
    assert approx(ratios[999], 5.0)  # fallback for uncalibrated size
    assert parse_bb_list(DEFAULT_BB) == [1280, 102400, 1048576]
    for bad in ("", ",", "1280,", "1280,,1024", "-1", "0", "10abc"):
        try:
            parse_bb_list(bad)
            raise SystemExit("self-test FAILED: accepted bad --bb %r" % bad)
        except ValueError:
            pass

    rows = load_rows([sim_path])
    assert len(rows) == 2
    dense, ldpc = rows
    # dense: base = 5*100 + 100 + 200 + 300 + 40(ge real) = 1140
    assert approx(dense["base"], 1140.0), dense["base"]
    assert dense["ge_src"] == "real"
    assert approx(dense["heavy"], 156.0)
    # total@ratio3 = 1140 + 3*156 = 1608; pessimistic doubles heavy: 2076
    assert approx(total_ops(dense, 3.0, False), 1608.0)
    assert approx(total_ops(dense, 3.0, True), 2076.0)
    # re-scored failure: pdf {0:0.9, 7:0.05, 13:0.05}
    assert approx(dense["fails"][6], 0.10)
    assert approx(dense["fails"][8], 0.05)
    assert approx(dense["fails"][12], 0.05)
    assert approx(dense["fails"][16], 0.00)
    # ldpcdense: base = 4*100 + 80 + 150 + 250 + 30 = 910; heavy 324
    assert approx(ldpc["base"], 910.0), ldpc["base"]
    assert approx(ldpc["fails"][12], 0.05) and approx(ldpc["fails"][16], 0.01)
    # crossover: (1140 - 910) / (324 - 156) = 230/168
    r, verdict = crossover_ratio(dense, ldpc, False)
    assert approx(r, 230.0 / 168.0), r
    assert verdict == "beats below this ratio"
    r2, _ = crossover_ratio(dense, ldpc, True)
    assert approx(r2, 230.0 / (2.0 * 168.0)), r2

    # Old layout: heavy derived as H*inact + H^2, ge falls back to estimate
    old = load_rows([old_path])[0]
    assert old["ge_src"] == "est"
    assert approx(old["base"], 5 * 100 + 100 + 200 + 300 + 50.0)
    assert approx(old["heavy"], 6 * 20.0 + 36.0)

    # New --max-inact layout: runaway rows have H-independent failure but no
    # meaningful cost means, so total_ops is deliberately unrankable.
    runaway = load_rows([runaway_path])[0]
    assert approx(runaway["runaway_rate"], 1.0)
    assert runaway["fails"][6] == 1.0 and runaway["fails"][16] == 1.0
    assert total_ops(runaway, 3.0, False) == float("inf")

    # H-independence assert must fire on an inconsistent fail_rate
    bad_path = os.path.join(tmpdir, "bad.csv")
    with open(bad_path, "w") as f:
        f.write(SELFTEST_SIM.replace("0.100000,1.000000,1.0,13",
                                     "0.300000,1.000000,1.0,13"))
    try:
        load_rows([bad_path])
        raise SystemExit("self-test FAILED: H-independence assert missing")
    except AssertionError:
        pass

    # End-to-end report generation must not throw
    import io
    buf = io.StringIO()
    print_report(rows, [1280, 102400], ratios, False, out=buf)
    text = buf.getvalue()
    assert "crossover K=100 oh=0" in text and "1.369" in text, text
    buf = io.StringIO()
    print_report(rows + [runaway], [1280], ratios, False, out=buf)
    text = buf.getvalue()
    assert "runaway" in text and "1.000000" in text, text

    print("rank_total self-test passed")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csvs", nargs="*", help="precode_sim result CSVs")
    ap.add_argument("--xor-csv", help="cold_xor_calibration.csv path")
    ap.add_argument("--muladd-csv", help="muladd_calibration.csv path")
    ap.add_argument("--bb", default=DEFAULT_BB,
                    help="comma list of block sizes (default %s)" % DEFAULT_BB)
    ap.add_argument("--pessimistic", action="store_true",
                    help="double the heavy GF(256) term")
    ap.add_argument("--assume-ratio", type=float, default=4.0,
                    help="muladd:xor ratio used when calibration is "
                         "missing for a block size (default 4.0)")
    ap.add_argument("--self-test", action="store_true",
                    help="run the built-in unit test against synthetic CSVs")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if not args.csvs:
        ap.error("no input CSVs (or use --self-test)")

    try:
        bb_list = parse_bb_list(args.bb)
    except ValueError as e:
        ap.error(str(e))
    xor_us = load_us_per_op(args.xor_csv, "xor")
    muladd_us = load_us_per_op(args.muladd_csv, "muladd")
    ratios = build_ratios(bb_list, xor_us, muladd_us, args.assume_ratio)
    rows = load_rows(args.csvs)
    if not rows:
        sys.stderr.write("rank_total: no data rows found\n")
        return 1
    print_report(rows, bb_list, ratios, args.pessimistic)
    return 0


if __name__ == "__main__":
    sys.exit(main())
