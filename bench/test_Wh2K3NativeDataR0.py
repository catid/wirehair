#!/usr/bin/env python3
"""Bounded native-fixture authentication tests; never invoke a codec/selector."""
import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('_k3_native_data_tests', HERE / 'Wh2K3NativeDataR0.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class DataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = M.load_report()

    def test_projection_and_fixed_lookup(self):
        data = M.extract(self.report)
        self.assertEqual(len(data['traces']), 6216)
        self.assertEqual(len(data['windows']), 43)
        self.assertEqual(len(data['rows']), 2226)
        self.assertEqual(sum(len(row['widths']) for row in data['history']), 53)
        self.assertEqual(sum(len(row['widths']) * len(row['ids']) for row in data['history']), 160)
        packed = M.build_lookup(data['pair'])
        self.assertEqual(len(packed), 13056)
        self.assertEqual(M.C.sha(packed), M.LOOKUP_SHA)
        header = M.render(data, packed)
        self.assertLess(len(header), M.CAP)
        self.assertIn(M.RAW_SHA.encode('ascii'), header)
        self.assertIn(b'kHistory[]', header)

    def test_tampered_member_rejected(self):
        original = M.C.read_regular

        def changed(path, *args, **kwargs):
            raw = original(path, *args, **kwargs)
            return raw + b' ' if path.name == 'raw.json' else raw

        with mock.patch.object(M.C, 'read_regular', side_effect=changed):
            with self.assertRaisesRegex(ValueError, 'member identity'):
                M.load_report()

    def test_projection_rejects_missing_and_changed_evidence(self):
        transforms = [lambda r: r['fresh'].pop(), lambda r: r['history'].pop(),
                      lambda r: r['hard'][0]['ranks'].__setitem__(0, 2),
                      lambda r: r['hard'][0]['ids'].__setitem__(0, r['hard'][0]['ids'][1]),
                      lambda r: r['pair'][0][0].__setitem__(2, 10),
                      lambda r: r['inputs']['prefixes'][0].__setitem__('original_widths', []),
                      lambda r: r['seams'][0].__setitem__('deficient', [[0, 1, 2]]),
                      lambda r: r['evidence']['unique_rows'].pop()]
        for change in transforms:
            report = copy.deepcopy(self.report)
            change(report)
            with self.assertRaises(ValueError):
                M.extract(report)

    def test_immutable_build_output(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k3-data-neutral-') as directory:
            path = Path(directory) / 'generated.inc'
            M.write_header(path, b'neutral')
            M.write_header(path, b'neutral')
            with self.assertRaisesRegex(ValueError, 'differs'):
                M.write_header(path, b'changed')
            self.assertEqual(path.read_bytes(), b'neutral')

    def test_bad_ids_and_pair(self):
        for values in ([0, 0, 2], [-1, 1, 2], [0, 1, 1 << 32], [False, 1, 2]):
            with self.assertRaises(ValueError): M.ids(values, 3)
        for pair in ([], [[[0] * 3] * 3], [[[0] * 3] * 3, [[0] * 2] * 3]):
            with self.assertRaises(ValueError): M.build_lookup(pair)


if __name__ == '__main__':
    unittest.main()
