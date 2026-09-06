#!/usr/bin/env python3
"""Neutral parser tests and read-only authentication; never score a candidate."""
import copy
import contextlib
import importlib.util
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest import mock


SPEC = importlib.util.spec_from_file_location(
    "_thue_recovery_history_under_test", Path(__file__).resolve().with_name(
        "Wh2ThueMorseRecoveryHistoryR0.py"))
H = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(H)


def origin(**changes):
    row = dict(issue=52, file="train.jsonl", line=2, K=6, B=2, arm="candidate", attempt=7,
               root="0x1234567890abcdef", schedule="repair-only", overhead=0, ids=list(range(6, 12)))
    row.update(changes)
    return row


def cell(**changes):
    row = dict(type="cell", K=6, B=2, arm="candidate", attempt=7,
               root="0x1234567890abcdef", schedule="repair-only", ids=list(range(6, 16)),
               overheads=[0, 1, 2, 4], outcomes=[1, 1, 0, 0])
    row.update(changes)
    return row


def jsonl(*rows):
    return b"".join(H.C.canonical(row) + b"\n" for row in rows)


class HistoryTests(unittest.TestCase):
    def bundle(self, structured=False):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        payload = b'{"neutral":true}\n'
        (directory / "data.json").write_bytes(payload)
        files = {"data.json": dict(bytes=len(payload), sha256=H.C.sha(payload)) if structured else H.C.sha(payload)}
        manifest = dict(protocol="wirehair.wh2.thue-morse-r0", outcome="PASS", files=files) if structured else files
        raw = H.C.canonical(manifest) + (b"" if structured else b"\n")
        (directory / "COMPLETE").write_bytes(raw)
        return directory, H.C.sha(raw), payload

    # Python 3.8 lacks TestCase.enterContext.
    def enterContext(self, context):
        result = context.__enter__()
        self.addCleanup(context.__exit__, None, None, None)
        return result

    def test_import_and_receipt_do_not_import_campaign_modules(self):
        self.assertNotIn("RADIX", vars(H))
        self.assertFalse(any(name in vars(H) for name in ("RowMapper", "fixed_feedback", "run_screen")))
        neutral = dict(provenance={"neutral": True}, prefixes=[{}], origins=[{}, {}], excluded_roots=[],
                       fixtures=[], pair=dict(pair_sha256="p", selection_sha256="s"))
        with mock.patch.object(H, "load_inputs", return_value=neutral) as load:
            receipt = H.input_receipt(123.0)
        load.assert_called_once_with(123.0)
        self.assertEqual(receipt["counts"], dict(prefixes=1, origins=2, excluded_roots=0, fixtures=0))
        self.assertEqual(receipt["inputs_sha256"], H.C.sha(H.C.canonical(neutral)))
        self.assertNotIn("A0", receipt)

    def test_both_manifest_formats_and_exact_roster(self):
        for structured in (False, True):
            directory, pin, payload = self.bundle(structured)
            budget = [0]
            content, receipt = H.read_bundle(directory, "COMPLETE", pin, ("data.json",), structured, budget)
            self.assertEqual(content, {"data.json": payload})
            self.assertEqual(receipt["files"]["data.json"], dict(bytes=len(payload), sha256=H.C.sha(payload)))
            self.assertEqual(budget[0], len(payload) + receipt["manifest_bytes"])
            (directory / "extra").write_bytes(b"")
            with self.assertRaises(ValueError):
                H.read_bundle(directory, "COMPLETE", pin, ("data.json",), structured, [0])

    def test_manifest_member_pin_canonical_and_caps(self):
        directory, pin, _ = self.bundle()
        with self.assertRaises(ValueError):
            H.read_bundle(directory, "COMPLETE", "0" * 64, ("data.json",), False, [0])
        with mock.patch.object(H, "INPUT_CAP", 1):
            with self.assertRaises(ValueError):
                H.read_bundle(directory, "COMPLETE", pin, ("data.json",), False, [0])
        with mock.patch.object(H, "AGGREGATE_CAP", 1):
            with self.assertRaises(ValueError):
                H.read_bundle(directory, "COMPLETE", pin, ("data.json",), False, [0])
        raw = (directory / "COMPLETE").read_bytes() + b"\n"
        (directory / "COMPLETE").write_bytes(raw)
        with self.assertRaises(ValueError):
            H.read_bundle(directory, "COMPLETE", H.C.sha(raw), ("data.json",), False, [0])
        raw = b'{"data.json":"x","data.json":"x"}\n'
        (directory / "COMPLETE").write_bytes(raw)
        with self.assertRaises(ValueError):
            H.read_bundle(directory, "COMPLETE", H.C.sha(raw), ("data.json",), False, [0])

    def test_directory_roster_stops_at_first_unexpected_entry(self):
        directory, pin, _ = self.bundle()

        def oversized_directory():
            yield Path("COMPLETE")
            yield Path("data.json")
            yield Path("unexpected")
            raise AssertionError("directory reader continued after unexpected member")

        with mock.patch.object(H.os, "scandir", return_value=contextlib.nullcontext(oversized_directory())):
            with self.assertRaisesRegex(ValueError, "directory roster"):
                H.read_bundle(directory, "COMPLETE", pin, ("data.json",), False, [0])

    def test_member_change_symlink_and_hardlink_rejected(self):
        directory, pin, _ = self.bundle()
        target = directory / "data.json"
        target.write_bytes(b"changed")
        with self.assertRaises(ValueError):
            H.read_bundle(directory, "COMPLETE", pin, ("data.json",), False, [0])
        elsewhere = Path(self.enterContext(tempfile.TemporaryDirectory())) / "target"
        elsewhere.write_bytes(b"neutral")
        target.unlink()
        target.symlink_to(elsewhere)
        with self.assertRaises(ValueError):
            H.read_bundle(directory, "COMPLETE", pin, ("data.json",), False, [0])
        target.unlink()
        os.link(str(elsewhere), str(target))
        with self.assertRaises(ValueError):
            H.read_bundle(directory, "COMPLETE", pin, ("data.json",), False, [0])

    def test_origin_order_first_witness_both_arms_and_prefix_union(self):
        witness = cell(overheads=[0], outcomes=[1], ids=list(range(6, 12)))
        attempt = dict(type="attempt", K=6, B=2, attempt=7, accepted=False, witness=witness)
        ignored = dict(type="attempt", accepted=False, witness=None)
        streams = {"validate.jsonl": jsonl(cell(arm="baseline", B=1280)),
                   "search.jsonl": jsonl(ignored, attempt),
                   "holdout.jsonl": jsonl(dict(type="begin"), cell())}
        result = H.extract_origins([origin(), dict(K=3)], streams)
        self.assertEqual([(r["file"], r["line"], r["overhead"]) for r in result],
                         [("train.jsonl", 2, 0), ("holdout.jsonl", 2, 0), ("holdout.jsonl", 2, 1),
                          ("search.jsonl", 2, 0), ("validate.jsonl", 1, 0), ("validate.jsonl", 1, 1)])
        prefixes = H.prefix_ledger(result)
        self.assertEqual(len(prefixes), 2)
        self.assertEqual([r["overhead"] for r in prefixes], [0, 1])
        self.assertEqual([r["original_widths"] for r in prefixes], [[2, 1280], [2, 1280]])
        self.assertEqual(H.prefix_ledger([origin(ids=list(reversed(range(6, 12))))])[0]["ids"],
                         list(reversed(range(6, 12))))

    def test_origin_types_coordinates_and_lengths_rejected(self):
        changes = ({"K": 6.0}, {"B": True}, {"attempt": False}, {"line": 0}, {"issue": 53},
                   {"file": "../x.jsonl"}, {"overhead": 1.0}, {"ids": [True] + list(range(7, 12))},
                   {"ids": [6] * 6}, {"ids": [1 << 32] + list(range(7, 12))},
                   {"root": "0x1234"}, {"root": None}, {"arm": "unknown"})
        for change in changes:
            with self.subTest(change=change), self.assertRaises(ValueError):
                H.validate_origin(origin(**change))

    def test_jsonl_bad_types_duplicate_keys_and_deadline(self):
        for row in (cell(K=6.0), cell(outcomes=[True, 1, 0, 0]), cell(overheads=[0, 1]),
                    cell(ids=[False] + list(range(7, 16)))):
            with self.assertRaises(ValueError):
                H.extract_origins([], {"holdout.jsonl": jsonl(row)})
        for raw in (b'{"type":"cell","type":"cell"}\n', b'{"x":NaN}\n', b'{}\n\n'):
            with self.assertRaises(ValueError):
                H.extract_origins([], {"holdout.jsonl": raw})
        with self.assertRaises(TimeoutError):
            H.extract_origins([], {"holdout.jsonl": jsonl(cell())}, time.monotonic() - 1)

    def test_recorded_fixture_extraction_no_rank_evaluation(self):
        final = [dict(family="neutral", index=0, ids=list(range(10)),
                      prefix_ranks={"6": 5, "7": 5, "8": 6, "10": 6}),
                 dict(family="neutral", index=1, ids=list(range(20, 30)),
                      prefix_ranks={"6": 5, "7": 6, "8": 6, "10": 6})]
        fixtures = H.extract_fixtures(final)
        self.assertEqual([(len(r["ids"]), r["expected_rank"]) for r in fixtures],
                         [(6, 5), (7, 5), (6, 5), (8, 6), (8, 6)])
        for field, value in (("6", 4), ("8", 5), ("6", 5.0), ("8", True)):
            malformed = copy.deepcopy(final)
            malformed[0]["prefix_ranks"][field] = value
            with self.assertRaises(ValueError):
                H.extract_fixtures(malformed)

    def test_actual_historical_bytes_authenticate_without_scoring(self):
        inputs = H.load_inputs(time.monotonic() + 8)
        self.assertEqual(set(inputs), {"provenance", "prefixes", "origins", "excluded_roots", "pair", "fixtures"})
        self.assertEqual([len(inputs[k]) for k in ("prefixes", "origins", "excluded_roots", "fixtures")],
                         [82, 187, 297, 33])
        self.assertEqual(H.C.sha(H.C.canonical(inputs["prefixes"])), H.PREFIX_SHA256)
        self.assertEqual(H.C.sha(H.C.canonical(inputs["origins"])), H.ORIGIN_SHA256)
        self.assertEqual(inputs["pair"]["pair_sha256"], H.PAIR_SHA256)
        self.assertEqual(len(set(H.MAIN_ROOTS)), 12)
        self.assertEqual(H.HARD_VALIDATION_ROOTS,
                         ("0xc0ac29b7c97c50dd", "0x3f84d5b5b5470917", "0x9216d5d98979fb1b"))
        self.assertTrue(set(H.MAIN_ROOTS) <= set(inputs["excluded_roots"]))
        self.assertEqual(sum(r["expected_rank"] == 5 for r in inputs["fixtures"]), 17)
        self.assertLess(len(H.C.canonical(inputs)), 100000)


if __name__ == "__main__":
    unittest.main()
