#!/usr/bin/env python3
"""Neutral .82 tests: no frozen pair selection or candidate recovery scoring."""
import contextlib
import hashlib
import importlib.util
import io
import itertools
from pathlib import Path
import struct
import sys
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = sibling("_k3_tm_test_worker", "Wh2K3ThueMorseR0.py")
W = sibling("_k3_tm_test_wrapper", "Wh2K3ThueMorseRunR0.py")
OLD = sibling("_k3_tm_test_old", "test_Wh2NoncommutingRadixRunR0.py")
OLD.M = W.C


class ArithmeticTests(unittest.TestCase):
    def test_all_field_products_and_determinants(self):
        M.F.init_field()
        for x, y in itertools.product(range(256), repeat=2):
            self.assertEqual(M.multiply_polynomial(x, y), M.F.MUL[x][y])
        for i in range(128):
            data = hashlib.sha256(("k3-neutral-matrix/" + str(i)).encode()).digest()
            matrix = [data[j:j + 3] for j in (0, 3, 6)]
            det = M.determinant3(matrix, M.multiply_polynomial)
            self.assertEqual(M.determinant3(matrix), det)
            self.assertEqual(det != 0, M.F.matrix_rank(matrix) == 3)
        self.assertEqual(M.checked_rank([(0, 0, 0)] * 3), 0)
        self.assertEqual(M.checked_rank([(1, 0, 0)] * 7), 1)
        self.assertEqual(M.checked_rank([(1, 0, 0), (0, 1, 0)] * 3), 2)

    def test_local_selector_stops_without_recovery(self):
        record, budget = [], mock.Mock()
        # Synthetic feedback, synthetic words and injected determinants only.
        with mock.patch.object(M, "determinant3", side_effect=[0, 1]), \
                mock.patch.object(M, "trace") as trace, \
                mock.patch.object(M, "history_inputs") as history:
            pair = M.choose_pair((3, 5, 11), ("0",), ((0, 1, 2),), budget, record)
        self.assertEqual([r["parameter"] for r in record], [1, 2])
        self.assertIsNotNone(record[0]["first_failure"])
        self.assertIsNone(record[1]["first_failure"])
        self.assertEqual(pair[1], M.companion((1, 5, 11)))
        trace.assert_not_called()
        history.assert_not_called()
        with mock.patch.object(M, "determinant3", return_value=0):
            record = []
            self.assertIsNone(M.choose_pair((3, 5, 11), ("0",), ((0, 1, 2),), budget, record))
        self.assertEqual(len(record), 254)
        self.assertNotIn(3, [r["parameter"] for r in record])

    def test_unrelated_lookup_pair_sequential_and_boundaries(self):
        for feedback in ((3, 5, 11), (3, 5, 11, 13)):
            pair = (M.companion(feedback), M.companion((19,) + feedback[1:]))
            mapper = M.Mapper(pair, M.F.Budget())
            n, product = len(feedback), M.identity(len(feedback))
            for packet_id in range(2051):
                self.assertEqual(mapper.row(packet_id), tuple(row[0] for row in product))
                product = M.F.matrix_multiply(product, pair[M.parity(packet_id)])
            for bit in range(2, 32):
                for offset in (-1, 0, 1):
                    packet_id = (1 << bit) + offset
                    self.assertEqual(mapper.row(packet_id), mapper.reference_row(packet_id))
            self.assertEqual(mapper.row(M.MAX_ID), mapper.reference_row(M.MAX_ID))
            self.assertEqual(tuple(mapper.row(i) for i in range(n)), M.identity(n))
            for invalid in (-1, 1 << 32, True, 1.0):
                with self.assertRaises(ValueError):
                    mapper.row(invalid)
            with mock.patch.object(mapper, "reference_row", return_value=(0,) * n):
                with self.assertRaisesRegex(ValueError, "disagreement"):
                    mapper.row(987654)

    def test_window_preserves_every_failure(self):
        mapper = mock.Mock()
        mapper.row.return_value = (1, 0, 0)
        row = M.check_window(range(7), mapper)
        self.assertEqual(row["deficient"], [list(t) for t in M.TRIPLES])
        self.assertEqual(mapper.row.call_count, 7)

    def test_main_contract_trace_manifest_without_equations(self):
        sha = hashlib.sha256()
        with mock.patch.object(M, "Mapper") as mapper, mock.patch.object(M, "choose_pair") as selector:
            for trial, root in enumerate(M.H.HARD_TRAINING_ROOTS):
                for schedule_index, schedule in enumerate(M.SCHEDULES):
                    ids = M.trace(2, root, schedule)
                    sha.update(struct.pack("<Q7I", trial * 120 + schedule_index * 30 + 1, *ids))
        self.assertEqual(sha.hexdigest(), M.DEVELOPMENT_TRACE_SHA)
        mapper.assert_not_called()
        selector.assert_not_called()
        for B, schedule in itertools.product(M.WIDTHS, M.SCHEDULES):
            ids = M.trace(B, "0x0000000000000000", schedule)
            self.assertEqual(len(set(ids)), 7)
            self.assertTrue(all(0 <= i <= M.MAX_ID for i in ids))
            if schedule == "adversarial":
                self.assertEqual(ids, sorted(ids, reverse=True))
            else:
                self.assertEqual(ids, sorted(ids))

    def test_invalid_cli_never_scores(self):
        for argv in ([], ["--worker", "extra"], ["--run"], ["--selftest"], ["--worker=1"]):
            with mock.patch.object(M, "run_screen") as run, contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(M.main(argv), 2)
                run.assert_not_called()

    def test_fresh_gate_every_cell_and_exact_boundary(self):
        rows = [dict(B=B, schedule=schedule, root=str(i), ranks=[3] * 5)
                for B, schedule in itertools.product(M.WIDTHS, M.SCHEDULES) for i in range(100)]
        self.assertTrue(M.summarize_fresh(rows, 100)["fresh_pass"])
        for offset in range(0, len(rows), 100):
            rows[offset]["ranks"] = [2, 3, 3, 3, 3]
            self.assertTrue(M.summarize_fresh(rows, 100)["fresh_pass"])
            rows[offset + 1]["ranks"] = [2, 2, 2, 2, 2]
            summary = M.summarize_fresh(rows, 100)
            self.assertFalse(summary["fresh_pass"])
            self.assertEqual(summary["cells"][offset // 100]["first_success"], [98, 1, 0, 0, 0, 1])
            rows[offset]["ranks"] = rows[offset + 1]["ranks"] = [3] * 5
        for corrupted in (rows[:-1], rows + rows[:1], rows[:-1] + rows[:1]):
            with self.assertRaises(ValueError):
                M.summarize_fresh(corrupted, 100)

    def test_structural_failure_never_enters_fresh_phase(self):
        # Same permutation matrix twice: unrelated synthetic cycle, not the
        # frozen feedback. All selection and scientific gates are injected.
        pair = (M.companion((1, 0, 0)),) * 2

        class FakeMapper:
            payload = bytes(13056)

            def __init__(self, *args):
                self.cache = {}

            def row(self, packet_id):
                value = M.identity(3)[packet_id % 3]
                self.cache[packet_id] = value
                return value

        def window(ids, mapper):
            return dict(ids=list(ids), deficient=[[0, 1, 2]])

        def trace_row(B, root, schedule, mapper):
            return dict(B=B, root=root, schedule=schedule, ids=list(range(7)), ranks=[3] * 5)

        with mock.patch.object(M, "history_inputs", return_value={"prefixes": [], "excluded_roots": []}), \
                mock.patch.object(M, "choose_pair", return_value=pair), \
                mock.patch.object(M, "Mapper", FakeMapper), \
                mock.patch.object(M, "checked_rank", return_value=3), \
                mock.patch.object(M, "check_window", side_effect=window), \
                mock.patch.object(M, "trace_result", side_effect=trace_row) as traces, \
                mock.patch.object(M, "summarize_fresh") as fresh:
            result = M.run_screen()
        self.assertEqual(result["outcome"], "FAIL", result.get("error"))
        self.assertFalse(result["summary"]["fresh_entered"])
        self.assertEqual(result["fresh"], [])
        self.assertEqual(traces.call_count, 72)
        fresh.assert_not_called()


class WrapperTests(unittest.TestCase):
    def test_inert_import_and_exact_identity(self):
        with mock.patch.object(W.C.subprocess, "Popen") as spawn:
            other = sibling("_k3_tm_reimport", "Wh2K3ThueMorseRunR0.py")
        spawn.assert_not_called()
        self.assertEqual(other.C.PROTOCOL, "wirehair.wh2.k3-thue-morse-r0")
        self.assertEqual(str(other.C.OUTPUT), "/var/tmp/wh2-k3-thue-morse-r0")
        self.assertEqual(other.C.SOURCES[0], "bench/Wh2K3ThueMorseR0.py")
        self.assertEqual(len(other.C.SOURCES), len(set(other.C.SOURCES)))
        self.assertIsNot(other.C, W.C)

    def test_receipt_pins_history_without_worker(self):
        with mock.patch.object(W, "_current_receipt", return_value={}), \
                mock.patch.object(W.H, "read_bundle", return_value=({}, {"authenticated": True})) as history, \
                mock.patch.object(M, "run_screen") as run:
            self.assertEqual(W.current_receipt(), {"history": {"authenticated": True}})
        history.assert_called_once()
        run.assert_not_called()


def load_tests(loader, tests, pattern):
    return unittest.TestSuite((tests, loader.loadTestsFromTestCase(OLD.FileTests),
                              loader.loadTestsFromTestCase(OLD.CaptureTests),
                              loader.loadTestsFromTestCase(OLD.PublicationTests)))


if __name__ == "__main__":
    unittest.main()
