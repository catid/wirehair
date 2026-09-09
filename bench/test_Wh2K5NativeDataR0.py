#!/usr/bin/env python3
"""Sealed K5 fixture projection tests; no codec, selector or fresh campaign."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

SPEC = importlib.util.spec_from_file_location('_k5_native_tests', Path(__file__).with_name('Wh2K5NativeDataR0.py'))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class DataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = M.load_report()
        cls.lookup = M.build_lookup()

    def test_exact_projection_and_lookup(self):
        data = M.extract(self.report)
        self.assertEqual(len(data['traces']), 6216)
        self.assertEqual(len(data['history']), 54)
        self.assertEqual(len(data['windows']), 30)
        self.assertEqual(len(data['rows']), 2270)
        self.assertEqual(len(self.lookup), 29440)
        self.assertEqual(M.C.sha(self.lookup), M.LOOKUP_SHA)
        self.assertEqual(self.lookup[:25], bytes(int(r == c) for r in range(5) for c in range(5)))
        origins = self.report['inputs']['origins']
        expected_widths = {}
        for origin in origins:
            expected_widths.setdefault(tuple(origin['ids']), set()).add(origin['b'])
        self.assertEqual({tuple(p['ids']): set(p['widths']) for p in data['history']}, expected_widths)
        self.assertEqual(sum(len(p['widths']) for p in data['history']), 57)
        self.assertEqual(sum(len(p['widths']) * len(p['ids']) for p in data['history']), 290)
        for actual, expected in zip(data['traces'], self.report['fresh'] + self.report['hard']):
            self.assertEqual(actual, {key: expected[key] for key in ('B', 'ids', 'ranks')})
        self.assertEqual(sum(t['ranks'][0] < 5 for t in data['traces']), 11)
        self.assertTrue(all(t['ranks'][1:] == [5] * 4 for t in data['traces']))

    def test_rendered_roster(self):
        header = M.render(M.extract(self.report), self.lookup)
        self.assertLess(len(header), M.CAP)
        self.assertIn(M.RAW_SHA.encode('ascii'), header)
        self.assertIn(M.LOOKUP_SHA.encode('ascii'), header)
        for name, count in ((b'Trace kTraces', 6216), (b'Prefix kHistory', 54),
                            (b'std::uint32_t kWindows', 30), (b'Row kRows', 2270)):
            body = header.split(b'static const ' + name, 1)[1].split(b' = {\n', 1)[1].split(b'};', 1)[0]
            self.assertEqual(len(body.splitlines()), count)

    def test_every_emitted_value_matches_sealed_evidence(self):
        header = M.render(M.extract(self.report), self.lookup).decode('ascii')

        def values(name):
            body = header.split(name, 1)[1].split(' = {\n', 1)[1].split('};', 1)[0]
            return [json.loads(line.rstrip(',').replace('{', '[').replace('}', ']').replace('u', ''))
                    for line in body.splitlines()]

        self.assertEqual(values('Trace kTraces[]'), [[r['B'], r['ids'], r['ranks']]
                         for r in self.report['fresh'] + self.report['hard']])
        self.assertEqual(values('Row kRows[]'), [[r['id'], r['row']]
                         for r in self.report['evidence']['unique_rows']])
        self.assertEqual(values('kWindows[][9]'), [r['ids'] for r in self.report['seams']])
        origins = self.report['inputs']['origins']
        expected = []
        for p in self.report['inputs']['prefixes']:
            widths = {r['b'] for r in origins if r['ids'] == p}
            expected.append([len(p), sum({2: 1, 64: 2, 1280: 4}[b] for b in widths), p])
        self.assertEqual(values('Prefix kHistory[]'), expected)
        body = header.split('kLookup[29440] = {\n', 1)[1].split('};', 1)[0]
        self.assertEqual(bytes(int(v) for v in body.replace('\n', '').rstrip(',').split(',')), self.lookup)

    def test_no_old_selection_or_receipt_revalidation(self):
        with mock.patch.object(M.R, 'load_report', side_effect=AssertionError('old bundle')):
            with mock.patch.object(M.C, 'current_receipt', side_effect=AssertionError('old HEAD')):
                self.assertEqual(M.extract(M.load_report()), M.extract(self.report))

    def test_tampered_bundle_rejected(self):
        original = M.C.read_regular
        for member in ('COMPLETE.json', 'raw.json', 'CLAIM.json', 'stderr.txt', 'summary.json'):
            def changed(path, *args, **kwargs):
                raw = original(path, *args, **kwargs)
                return raw + b' ' if path.name == member else raw
            with self.subTest(member=member), mock.patch.object(M.C, 'read_regular', side_effect=changed):
                with self.assertRaisesRegex(ValueError, 'identity'): M.load_report()

    def test_projection_rejects_bad_evidence(self):
        mutations = [lambda r: r['fresh'].pop(), lambda r: r['hard'].pop(),
                     lambda r: r['history'].pop(), lambda r: r['inputs']['origins'].pop(),
                     lambda r: r['fresh'][0]['ranks'].__setitem__(0, 4),
                     lambda r: r['hard'][0]['ranks'].__setitem__(0, 4),
                     lambda r: r['hard'][0]['ranks'].__setitem__(1, 4),
                     lambda r: r['hard'][0]['ids'].__setitem__(0, r['hard'][0]['ids'][1]),
                     lambda r: r['hard'][0].__setitem__('B', True),
                     lambda r: r['inputs']['origins'][0].__setitem__('b', 3),
                     lambda r: r['inputs']['prefixes'][0].__setitem__(0, 100000),
                     lambda r: r['history'][0].__setitem__('rank', 4),
                     lambda r: r['pair'][0][0].__setitem__(4, 120),
                     lambda r: r['seams'][0].__setitem__('deficient', [[0,1,2,3,4]]),
                     lambda r: r['seams'].pop(), lambda r: r['evidence']['unique_rows'].pop(),
                     lambda r: r['evidence']['unique_rows'][0]['row'].__setitem__(0, False)]
        for index, change in enumerate(mutations):
            r = copy.deepcopy(self.report)
            change(r)
            with self.subTest(index=index), self.assertRaises(ValueError): M.extract(r)

    def test_bad_lookup_rejected(self):
        for lookup in (self.lookup[:-1], bytes([self.lookup[0] ^ 1]) + self.lookup[1:]):
            with self.assertRaisesRegex(ValueError, 'lookup bytes/hash'):
                M.render(M.extract(self.report), lookup)

    def test_immutable_output_and_symlink_rejection(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k5-data-neutral-') as directory:
            path = Path(directory) / 'fixture.inc'
            M.R.write_header(path, b'neutral')
            M.R.write_header(path, b'neutral')
            with self.assertRaisesRegex(ValueError, 'differs'): M.R.write_header(path, b'changed')
            link = Path(directory) / 'link.inc'
            link.symlink_to(path)
            with self.assertRaisesRegex(ValueError, 'regular file'): M.R.write_header(link, b'neutral')
            self.assertEqual(path.read_bytes(), b'neutral')


if __name__ == '__main__':
    unittest.main()
