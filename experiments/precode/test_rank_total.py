#!/usr/bin/env python3
"""Unit, property, oracle, and repository-data tests for rank_total.py."""

import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import random
import statistics
import subprocess
import sys
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "rank_total.py"
SPEC = importlib.util.spec_from_file_location("rank_total", SCRIPT)
rank_total = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rank_total)


def base_row(scheme="dense", heavy_rows=6, trials=100,
             pdf="0:0.900000|7:0.100000"):
    bins = {int(item.split(":")[0]): float(item.split(":")[1])
            for item in pdf.split("|")}
    fail_rate = sum(p for d, p in bins.items() if d > heavy_rows)
    noheavy = sum(p for d, p in bins.items() if d > 0)
    completed = [d for d in bins if d != 999999]
    completed_mass = sum(p for d, p in bins.items() if d != 999999)
    def_mu = (sum(d * p for d, p in bins.items() if d != 999999)
              / completed_mass if completed_mass else 0.0)
    return {
        "K": "100", "scheme": scheme, "D": "4", "H": str(heavy_rows),
        "oh": "0", "trials": str(trials),
        "fail_rate": "%.6f" % fail_rate,
        "fail_rate_noheavy": "%.6f" % noheavy,
        "def_mu": "%.6f" % def_mu,
        "def_max": str(max((d for d in completed if bins[d] > 0), default=0)),
        "def_pdf": pdf, "inact_mu": "20.000000", "inact_sd": "1.0",
        "inact_max": "25", "rank_mu": "19.0", "peeled_mu": "80.0",
        "residual_rows_mu": "19.0", "recv_xors_per_packet": "5.0",
        "precode_gen_xors_mu": "100.0", "sparse_solve_xors_mu": "200.0",
        "backsub_xors_mu": "300.0", "ge_block_xors_mu": "50.0",
        "ge_bitops_mu": "999.0",
    }


def write_result(path, rows, header=rank_total.BASE_FIELDS):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def xor_calibration_row(block_bytes="1280", timing="1.0"):
    return [block_bytes, "1", "1", "1", "1.0", timing, "1.0", "0x1"]


def write_csv(path, header, rows, preamble=None):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        if preamble:
            handle.write(preamble + "\n")
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def default_options(**overrides):
    values = {
        "confidence": 0.95,
        "cost_h": "modeled",
        "max_observed_fail": 1.0,
        "max_upper_fail": 1.0,
        "min_trials": 1,
        "pareto_metric": "wilson_upper",
        "pessimistic": False,
        "reliability_h": "modeled",
    }
    values.update(overrides)
    return values


class PdfAndNumericTests(unittest.TestCase):
    def test_property_normalized_empirical_pdfs(self):
        rng = random.Random(0x52414E4B)
        for _ in range(250):
            trials = 10000
            cuts = sorted(rng.sample(range(1, trials), 4))
            counts = [cuts[0], cuts[1] - cuts[0], cuts[2] - cuts[1],
                      cuts[3] - cuts[2], trials - cuts[3]]
            deficits = rng.sample(range(0, 40), len(counts))
            text = "|".join("%d:%.6f" % (d, count / trials)
                            for d, count in zip(deficits, counts))
            pdf, recovered, _ = rank_total.parse_def_pdf(text, trials=trials)
            self.assertEqual(sum(recovered.values()), trials)
            self.assertAlmostEqual(math.fsum(pdf.values()), 1.0)
            for deficit, count in zip(deficits, counts):
                self.assertEqual(recovered[deficit], count)

    def test_bad_pdf_corpus(self):
        bad = [
            "", "0:nan|1:1", "0:inf|1:0", "-1:1", "65536:1",
            "0:0.5|0:0.5", "0:0.4|1:0.4", "0:1.1", "0:-0.1|1:1.1",
            "0", "0:1|", ":1", "0:", "0:0.333333|1:0.666667",
        ]
        for text in bad:
            with self.subTest(text=text), self.assertRaises(rank_total.RankError):
                trials = 10 if text == "0:0.333333|1:0.666667" else None
                rank_total.parse_def_pdf(text, "bad.csv", 2, trials)

    def test_parser_overflow_and_cli_numbers(self):
        huge = "9" * 10000
        with self.assertRaises(rank_total.RankError):
            rank_total.parse_uint(huge, "x", 1, "n")
        for text in ("nan", "inf", "-inf", "-1", "0"):
            with self.subTest(text=text), self.assertRaises(rank_total.RankError):
                rank_total.parse_cli_float(text, "ratio", 0.0, positive=True)

    def test_wilson_against_independent_oracle(self):
        # This transcription deliberately does not call the production helper.
        confidence = 0.95
        z = statistics.NormalDist().inv_cdf(confidence)
        for trials in (1, 2, 10, 100, 1000):
            for failures in sorted({0, trials // 3, trials}):
                p = failures / trials
                denominator = trials + z * z
                oracle = ((trials * p + z * z / 2
                           + z * math.sqrt(trials * p * (1 - p) + z * z / 4))
                          / denominator)
                actual = rank_total.wilson_upper(failures, trials, confidence)
                self.assertAlmostEqual(actual, oracle, places=14)
        self.assertAlmostEqual(rank_total.wilson_upper(0, 100, 0.95),
                               0.026342720783174303)


class ParserAndCostTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_exact_result_schema_and_cost(self):
        path = self.directory / "valid.csv"
        write_result(path, [base_row()])
        row = rank_total.load_rows([str(path)])[0]
        self.assertEqual(row["ge_source"], "estimated")
        self.assertEqual(row["heavy_source"], "derived")
        self.assertAlmostEqual(row["base_ops"], 1150.0)
        self.assertAlmostEqual(row["heavy_muladds"], 156.0)

        modeled = rank_total.evaluate_row(row, 1280, 3.0, default_options())
        self.assertAlmostEqual(modeled["total_ops"], 1618.0)
        self.assertEqual(modeled["cost_semantics"], "modeled_h_ledger")

        alternate = rank_total.evaluate_row(
            row, 1280, 3.0,
            default_options(reliability_h=12, cost_h="reliability"))
        self.assertEqual(alternate["failures"], 0)
        self.assertEqual(alternate["cost_h"], 12)
        self.assertAlmostEqual(alternate["heavy_muladds"], 384.0)
        self.assertAlmostEqual(alternate["total_ops"], 2302.0)
        self.assertEqual(alternate["heavy_source"], "recomputed")
        self.assertIn("base_held", alternate["cost_semantics"])

    def test_property_cost_calculations(self):
        path = self.directory / "property.csv"
        write_result(path, [base_row(pdf="0:1.000000")])
        template = rank_total.load_rows([str(path)])[0]
        rng = random.Random(0x434F5354)
        for _ in range(250):
            row = dict(template)
            row["base_ops"] = rng.uniform(0.0, 1e8)
            row["heavy_muladds"] = rng.uniform(0.0, 1e5)
            row["inact_mu"] = rng.uniform(0.0, 1e4)
            ratio = rng.uniform(0.001, 100.0)
            pessimistic = bool(rng.randrange(2))
            factor = 2.0 if pessimistic else 1.0

            modeled = rank_total.evaluate_row(
                row, 1280, ratio,
                default_options(pessimistic=pessimistic))
            self.assertAlmostEqual(
                modeled["total_ops"],
                row["base_ops"] + ratio * factor * row["heavy_muladds"])

            alternate_h = rng.randrange(129)
            recomputed = rank_total.evaluate_row(
                row, 1280, ratio,
                default_options(pessimistic=pessimistic,
                                reliability_h=alternate_h,
                                cost_h="reliability"))
            oracle_heavy = alternate_h * row["inact_mu"] + alternate_h ** 2
            self.assertAlmostEqual(
                recomputed["total_ops"],
                row["base_ops"] + ratio * factor * oracle_heavy)

    def test_recorded_zero_costs_do_not_trigger_fallbacks(self):
        header = rank_total.BASE_FIELDS + rank_total.HEAVY_FIELDS \
            + rank_total.REPLAY_FIELDS
        row = base_row()
        row.update({name: "0" for name in rank_total.HEAVY_FIELDS
                    + rank_total.REPLAY_FIELDS})
        row["heavy_divs_mu"] = row["H"]
        path = self.directory / "zero.csv"
        write_result(path, [row], header)
        parsed = rank_total.load_rows([str(path)])[0]
        self.assertEqual(parsed["ge_source"], "real")
        self.assertEqual(parsed["heavy_source"], "recorded")
        self.assertEqual(parsed["heavy_muladds"], 0.0)
        self.assertAlmostEqual(parsed["base_ops"], 1100.0)

    def test_modeled_h_crossover_is_preserved(self):
        dense_path = self.directory / "dense.csv"
        alternate_path = self.directory / "alternate.csv"
        dense = base_row("dense", heavy_rows=6, pdf="0:1.000000")
        alternate = base_row("ldpcdense_s8_d4_h12", heavy_rows=12,
                             pdf="0:1.000000")
        alternate["precode_gen_xors_mu"] = "0.0"
        write_result(dense_path, [dense])
        write_result(alternate_path, [alternate])
        rows = rank_total.load_rows([str(dense_path), str(alternate_path)])
        dense_row = next(row for row in rows if row["scheme"] == "dense")
        alternate_row = next(row for row in rows
                             if row["scheme"].startswith("ldpcdense"))
        ratio, verdict = rank_total.crossover_ratio(dense_row, alternate_row)
        self.assertAlmostEqual(ratio, 100.0 / (384.0 - 156.0))
        self.assertEqual(verdict, "beats below this ratio")

    def test_runaway_cost_is_unavailable_and_never_ranked(self):
        header = (rank_total.BASE_FIELDS + rank_total.HEAVY_FIELDS
                  + ("runaway_rate",))
        row = base_row(pdf="999999:1.000000")
        row.update({"heavy_muladds_mu": "0", "heavy_divs_mu": "0",
                    "runaway_rate": "1.000000"})
        path = self.directory / "runaway.csv"
        write_result(path, [row], header)
        parsed = rank_total.load_rows([str(path)])[0]
        evaluated = rank_total.evaluate_row(
            parsed, 1280, 1.0, default_options())
        self.assertFalse(evaluated["eligible"])
        self.assertIsNone(evaluated["total_ops"])
        self.assertIn("runaway_rate>0", evaluated["rejection_reasons"])

        mismatch = base_row(trials=100000,
                            pdf="0:0.999990|999999:0.000010")
        mismatch.update({"heavy_muladds_mu": "0", "heavy_divs_mu": "0",
                         "runaway_rate": "0.000000"})
        mismatch_path = self.directory / "runaway_mismatch.csv"
        write_result(mismatch_path, [mismatch], header)
        with self.assertRaisesRegex(rank_total.RankError,
                                    r"runaway_mismatch.csv:2: runaway_rate"):
            rank_total.load_rows([str(mismatch_path)])

    def test_strict_calibration_schemas_and_preamble(self):
        path = self.directory / "xor.csv"
        write_csv(path, rank_total.XOR_CALIBRATION_FIELDS,
                  [xor_calibration_row()], "# command=test")
        self.assertEqual(rank_total.load_calibration(str(path), "xor"),
                         {1280: 1.0})

        mul = self.directory / "mul.csv"
        values = ["1280", "1", "1", "1", "1.0", "1.0", "2.0",
                  "1.0", "2.0000", "0x2"]
        write_csv(mul, rank_total.MULADD_CALIBRATION_FIELDS, [values])
        self.assertEqual(rank_total.load_calibration(str(mul), "muladd"),
                         {1280: 2.0})
        self.assertEqual(rank_total.load_calibration(str(mul), "xor"),
                         {1280: 1.0})

    def test_bad_calibration_corpus(self):
        mutations = [
            (5, "nan"), (5, "inf"), (5, "0"), (1, "-1"),
            (0, str(1 << 40)), (7, "bad-checksum"),
        ]
        for index, value in mutations:
            with self.subTest(index=index, value=value):
                path = self.directory / ("bad_%d.csv" % index)
                row = xor_calibration_row()
                row[index] = value
                write_csv(path, rank_total.XOR_CALIBRATION_FIELDS, [row])
                with self.assertRaises(rank_total.RankError):
                    rank_total.load_calibration(str(path), "xor")

        duplicate = self.directory / "duplicate.csv"
        write_csv(duplicate, rank_total.XOR_CALIBRATION_FIELDS,
                  [xor_calibration_row(), xor_calibration_row()])
        with self.assertRaises(rank_total.RankError):
            rank_total.load_calibration(str(duplicate), "xor")

        muladd_row = ["1280", "1", "1", "1", "1.0", "1.0", "2.0",
                      "1.0", "2.0000", "0x2"]
        for label, header, original in (
                ("xor", rank_total.XOR_CALIBRATION_FIELDS,
                 xor_calibration_row()),
                ("muladd", rank_total.MULADD_CALIBRATION_FIELDS,
                 muladd_row)):
            for index, field in enumerate(header):
                if field == "checksum":
                    continue
                for invalid in ("nan", "-1"):
                    with self.subTest(label=label, field=field,
                                      invalid=invalid):
                        row = list(original)
                        row[index] = invalid
                        path = self.directory / (
                            "%s_%s_%s.csv" % (label, field, invalid))
                        write_csv(path, header, [row])
                        with self.assertRaises(rank_total.RankError):
                            rank_total.load_calibration(str(path), label)

    def test_missing_extra_fields_and_duplicate_rows(self):
        missing = self.directory / "missing.csv"
        write_csv(missing, rank_total.BASE_FIELDS, [["1"]])
        with self.assertRaisesRegex(rank_total.RankError, r"missing.csv:2"):
            rank_total.load_rows([str(missing)])

        extra = self.directory / "extra.csv"
        write_result(extra, [base_row()])
        with open(extra, "a", encoding="utf-8") as handle:
            handle.write(",extra\n")
        with self.assertRaisesRegex(rank_total.RankError, r"extra.csv:3"):
            rank_total.load_rows([str(extra)])

        duplicate = self.directory / "duplicate.csv"
        write_result(duplicate, [base_row(), base_row()])
        with self.assertRaisesRegex(rank_total.RankError,
                                    r"duplicate.csv:3: duplicate"):
            rank_total.load_rows([str(duplicate)])

    def test_every_numeric_result_field_rejects_nonfinite_and_negative(self):
        header = (("pivot_window",) + rank_total.BASE_FIELDS
                  + rank_total.HEAVY_FIELDS + rank_total.REPLAY_FIELDS
                  + ("runaway_rate",))
        valid = base_row()
        valid.update({
            "pivot_window": "0", "heavy_muladds_mu": "156.0",
            "heavy_divs_mu": "6.0", "ge_real_bitops_mu": "400.0",
            "ge_real_rowops_mu": "40.0", "fill_in_mu": "2.0",
            "def_outside_w18_rate": "0.0", "def_band_w95": "3",
            "def_band_w99": "5", "runaway_rate": "0.0",
        })
        for field in header:
            if field in {"scheme", "def_pdf"}:
                continue
            for suffix, invalid in (("nan", "nan"), ("negative", "-1")):
                with self.subTest(field=field, invalid=invalid):
                    row = dict(valid)
                    row[field] = invalid
                    path = self.directory / ("%s_%s.csv" % (field, suffix))
                    write_result(path, [row], header)
                    with self.assertRaises(rank_total.RankError):
                        rank_total.load_rows([str(path)])


class ReliabilityAndParetoTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_pareto_against_pairwise_oracle_and_ties(self):
        rng = random.Random(321)
        for count in range(1, 30):
            rows = []
            for index in range(count):
                rows.append({
                    "row_id": str(index),
                    "total_ops": float(rng.randrange(5)),
                    "upper_failure": float(rng.randrange(5)),
                })
            expected = []
            for candidate in rows:
                if not any(
                    other is not candidate
                    and other["total_ops"] <= candidate["total_ops"]
                    and other["upper_failure"] <= candidate["upper_failure"]
                    and (other["total_ops"] < candidate["total_ops"]
                         or other["upper_failure"] < candidate["upper_failure"])
                    for other in rows
                ):
                    expected.append(candidate["row_id"])
            self.assertEqual(rank_total.pareto_ids(rows), sorted(expected))

        tied = [
            {"row_id": "a", "total_ops": 1.0, "upper_failure": 0.1},
            {"row_id": "b", "total_ops": 1.0, "upper_failure": 0.1},
        ]
        self.assertEqual(rank_total.pareto_ids(tied), ["a", "b"])

    def test_gates_small_zero_sample_insufficient_and_no_feasible(self):
        small_path = self.directory / "small.csv"
        write_result(small_path, [base_row(trials=10, pdf="0:1.000000")])
        row = rank_total.load_rows([str(small_path)])[0]
        evaluated = rank_total.evaluate_row(
            row, 1280, 1.0,
            default_options(min_trials=100, max_observed_fail=0.0,
                            max_upper_fail=0.05))
        self.assertFalse(evaluated["eligible"])
        self.assertIn("trials<100", evaluated["rejection_reasons"])
        self.assertTrue(any(reason.startswith("wilson_upper>")
                            for reason in evaluated["rejection_reasons"]))

        report = rank_total.build_report(
            [row], [1280], {1280: {"value": 1.0, "source": "assumed"}},
            default_options(min_trials=100), [str(small_path)])
        group = report["groups"][0]
        self.assertEqual(group["winner_ids"], [])
        self.assertEqual(group["pareto_ids"], [])
        self.assertFalse(group["rows"][0]["winner"])

    def test_winners_never_include_rejected_rows_and_cost_ties_co_win(self):
        path = self.directory / "rows.csv"
        good1 = base_row("good1", trials=100, pdf="0:1.000000")
        good2 = base_row("good2", trials=100, pdf="0:1.000000")
        bad = base_row("bad", trials=100, pdf="0:0.900000|7:0.100000")
        write_result(path, [good1, good2, bad])
        rows = rank_total.load_rows([str(path)])
        report = rank_total.build_report(
            rows, [1280], {1280: {"value": 1.0, "source": "assumed"}},
            default_options(max_observed_fail=0.01), [str(path)])
        group = report["groups"][0]
        winners = [row for row in group["rows"] if row["winner"]]
        self.assertEqual({row["scheme"] for row in winners}, {"good1", "good2"})
        self.assertTrue(all(row["eligible"] for row in winners))
        self.assertEqual({row["eligible_rank"] for row in winners}, {1})

        # Equal cost with worse reliability is dominated, not a co-winner.
        relaxed = rank_total.build_report(
            rows, [1280], {1280: {"value": 1.0, "source": "assumed"}},
            default_options(), [str(path)])
        relaxed_rows = relaxed["groups"][0]["rows"]
        bad_result = next(row for row in relaxed_rows if row["scheme"] == "bad")
        self.assertFalse(bad_result["winner"])
        self.assertFalse(bad_result["pareto"])


class CliAndRepositoryE2ETests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                              text=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, check=False)

    def test_invalid_inputs_emit_no_ranking_or_json(self):
        base = base_row()
        mutations = [
            ("fail_rate", "nan"), ("fail_rate", "inf"),
            ("fail_rate", "-1"), ("fail_rate", "1.1"),
            ("precode_gen_xors_mu", "-1"), ("inact_mu", "nan"),
            ("def_mu", "9"), ("def_max", "8"),
            ("inact_max", "1"), ("rank_mu", "99"),
            ("K", "4294967296"), ("trials", str(rank_total.MAX_TRIALS + 1)),
            ("def_pdf", "0:0.5|0:0.5"), ("def_pdf", "-1:1"),
            ("def_pdf", "0:0.9"), ("def_pdf", "0:0.5|7:nan"),
        ]
        for index, (field, value) in enumerate(mutations):
            with self.subTest(field=field, value=value):
                row = dict(base)
                row[field] = value
                path = self.directory / ("invalid_%02d.csv" % index)
                output = self.directory / ("invalid_%02d.json" % index)
                write_result(path, [row])
                output.write_text("sentinel", encoding="utf-8")
                result = self.run_cli(path, "--bb", "1280", "--json-out", output)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(result.stdout, "")
                self.assertEqual(output.read_text(encoding="utf-8"), "sentinel")
                self.assertRegex(result.stderr, r"invalid_\d+\.csv:2")

        text_report = self.directory / "generated-report.txt"
        text_report.write_text("rank_total report v1\n", encoding="utf-8")
        rejected_report = self.run_cli(text_report, "--bb", "1280")
        self.assertNotEqual(rejected_report.returncode, 0)
        self.assertEqual(rejected_report.stdout, "")
        self.assertIn("unsupported precode result header", rejected_report.stderr)

        valid = self.directory / "valid.csv"
        write_result(valid, [base_row()])

        malformed = self.directory / "wrong-fields.csv"
        write_csv(malformed, rank_total.BASE_FIELDS, [["1"]])
        malformed_run = self.run_cli(malformed, "--bb", "1280")
        self.assertNotEqual(malformed_run.returncode, 0)
        self.assertEqual(malformed_run.stdout, "")
        self.assertIn("wrong field count", malformed_run.stderr)

        bad_calibration = self.directory / "bad-calibration.csv"
        calibration_row = xor_calibration_row()
        calibration_row[5] = "nan"
        write_csv(bad_calibration, rank_total.XOR_CALIBRATION_FIELDS,
                  [calibration_row])
        calibration_output = self.directory / "bad-calibration.json"
        calibration_output.write_text("sentinel", encoding="utf-8")
        calibration_run = self.run_cli(
            valid, "--bb", "1280", "--xor-csv", bad_calibration,
            "--json-out", calibration_output)
        self.assertNotEqual(calibration_run.returncode, 0)
        self.assertEqual(calibration_run.stdout, "")
        self.assertEqual(calibration_output.read_text(encoding="utf-8"),
                         "sentinel")
        self.assertRegex(calibration_run.stderr,
                         r"bad-calibration.csv:2:.*finite")

        for option, value in (("--assume-ratio", "nan"),
                              ("--assume-ratio", "0"),
                              ("--confidence", "1"),
                              ("--max-upper-fail", "inf"),
                              ("--reliability-h", "129")):
            with self.subTest(option=option):
                result = self.run_cli(valid, "--bb", "1280", option, value)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(result.stdout, "")

    def test_deterministic_human_and_json_reports(self):
        first = self.directory / "a.csv"
        second = self.directory / "b.csv"
        write_result(first, [base_row("a", trials=100, pdf="0:1.000000")])
        write_result(second, [base_row("b", trials=100, pdf="0:1.000000")])
        out1 = self.directory / "one.json"
        out2 = self.directory / "two.json"
        args = ["--bb", "1280", "--min-trials", "100",
                "--max-upper-fail", "0.05"]
        run1 = self.run_cli(first, second, *args, "--json-out", out1)
        run2 = self.run_cli(second, first, *args, "--json-out", out2)
        self.assertEqual(run1.returncode, 0, run1.stderr)
        self.assertEqual(run2.returncode, 0, run2.stderr)
        self.assertEqual(run1.stdout, run2.stdout)
        self.assertEqual(out1.read_bytes(), out2.read_bytes())
        machine = json.loads(out1.read_text(encoding="utf-8"))
        self.assertEqual(machine["schema"], "wirehair.precode-ranking.v1")
        self.assertEqual(len(machine["groups"][0]["winner_ids"]), 2)

    def test_cli_defaults_reject_unreliable_cost_winner(self):
        path = self.directory / "unreliable.csv"
        output = self.directory / "unreliable.json"
        write_result(path, [base_row()])
        run = self.run_cli(path, "--bb", "1280", "--json-out", output)
        self.assertEqual(run.returncode, 0, run.stderr)
        report = json.loads(output.read_text(encoding="utf-8"))
        group = report["groups"][0]
        self.assertEqual(group["winner_ids"], [])
        self.assertFalse(group["rows"][0]["eligible"])
        self.assertIn("winner: no feasible scheme", run.stdout)

    def test_degree_rowdist_repository_conclusion(self):
        result_dir = HERE / "results"
        names = [
            "degree_rowdist_wirehair_K10000_32000_t500_20260708.csv",
            "degree_rowdist_wirehair_K64000_t100_20260708.csv",
            "degree_rowdist_lt_m2_c1024_K10000_32000_t500_20260708.csv",
            "degree_rowdist_lt_m2_c1024_K64000_t100_20260708.csv",
            "degree_rowdist_lt_m1_c16_K10000_32000_t500_20260708.csv",
            "degree_rowdist_lt_m1_c64_K10000_32000_t500_20260708.csv",
            "degree_rowdist_rs_c001_d50_c128_K10000_32000_t500_20260708.csv",
        ]
        paths = [str(result_dir / name) for name in names]
        rows = rank_total.load_rows(paths)
        xor_us = rank_total.load_calibration(
            str(HERE.parent / "peeling/results/cold_xor_calibration.csv"), "xor")
        muladd_us = rank_total.load_calibration(
            str(HERE.parent / "peeling/results/muladd_calibration.csv"), "muladd")
        ratios = rank_total.build_ratios([1280], xor_us, muladd_us, 4.0)
        report = rank_total.build_report(
            rows, [1280], ratios,
            default_options(min_trials=100, max_observed_fail=0.005,
                            max_upper_fail=0.05), paths)

        for group in report["groups"]:
            for row in group["rows"]:
                basename = os.path.basename(row["source"])
                if any(token in basename for token in
                       ("lt_m1_c16", "lt_m1_c64", "rs_c001_d50_c128")):
                    self.assertFalse(row["eligible"], basename)

        def codecport_rows(k, oh):
            group = next(group for group in report["groups"]
                         if (group["K"], group["oh"], group["block_bytes"])
                         == (k, oh, 1280))
            return [row for row in group["rows"]
                    if row["eligible"] and row["scheme"].startswith("codecport")]

        for k in (32000, 64000):
            candidates = codecport_rows(k, 0)
            production = [row for row in candidates
                          if "rowdist_wirehair" in row["source"]]
            retuned = [row for row in candidates
                       if "rowdist_lt_m2_c1024" in row["source"]]
            self.assertTrue(production and retuned)
            self.assertLess(min(row["total_ops"] for row in production),
                            min(row["total_ops"] for row in retuned))


if __name__ == "__main__":
    unittest.main()
