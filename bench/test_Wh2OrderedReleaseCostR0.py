"""Synthetic timing/ledger tests only: never invokes a scientific codec worker."""
import copy
import json
from pathlib import Path
import tempfile
import time
import unittest
import subprocess
from unittest.mock import patch

import Wh2OrderedReleaseCostR0 as C


def synthetic():
    fixtures, old = [], []
    for b in C.WIDTHS:
        arm = dict(profile='00'*32, packets='00'*(18*b), steps=[6, 6])
        fixtures.append(dict(width=b, source=bytes((37*i+i//11) % 256 for i in range(6*b)).hex(),
                             arms=[copy.deepcopy(arm) for _ in range(3)]))
        old.append(dict(handles=[dict(arm=3, profile_hex=arm['profile'])], packets_hex=[arm['packets']]*4))
    prelude = dict(clocks=[100, 100, 102, 112, 108, 120], before=[0]*4, after=[0]*4)
    header = dict(type='header', protocol=C.PROTOCOL, claim='0'*64, batch=C.BATCH,
                  identity_hex='00', prelude=prelude, fixtures=fixtures)
    records = []; previous = prelude; total = 0
    for coordinate in C.roster():
        metric, arm, q = coordinate[4], coordinate[7], coordinate[8]
        ready, cpu = previous['clocks'][5]+100, previous['clocks'][4]
        target = ready+q; duration = (100000, 99000, 120000)[arm]
        observation = dict(clocks=[target+4, cpu+40, target+6, target+6+duration,
                                  cpu+40+duration, target+20+duration], before=[0]*4, after=[0]*4)
        counts = [0 if metric else 128, 0 if metric else 2304, 128 if metric else 0,
                  768 if metric else 0, 128 if metric else 0, 128]
        records.append(dict(type='record', coordinate=coordinate, ready=ready, target=target,
                            wait=[ready+10, cpu+10, target+3, cpu+30], observation=observation,
                            counts=counts, addresses=[4096]*128, address_count=128, complete=True, checked=True))
        total += duration; previous = observation
    footer = dict(type='footer', complete=True, records=C.CALLBACKS, work_ns=total)
    return header, records, footer, dict(fixtures=old, identity_before=dict(canonical_hex='00'))


class CostTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.header, cls.records, cls.footer, cls.old = synthetic()

    def test_whole_synthetic_stream(self):
        raw = b''.join(C.A.canonical(row) for row in [self.header]+self.records+[self.footer])
        result = C.verify(raw, '0'*64, self.old)
        self.assertEqual(result['outcome'], 'PASS')
        self.assertTrue(result['all_WH1_faster'])
        self.assertEqual(len(result['statistics']), 90)
        self.assertLess(len(raw), C.RAW_CAP)

    def test_roster_independent_coverage(self):
        from collections import Counter
        counts, phases = Counter(), Counter()
        for i, c in enumerate(C.roster()):
            index, rep, order, width, metric, comparison, p, arm, q = c
            self.assertEqual(i, index)
            self.assertEqual(arm, C.PAIRS[comparison][int('010110100110010110'[p]) ^ order])
            if p == 0:
                counts[width, metric, comparison, order, rep] += 1
            if p >= 2 and p % 2 == 0:
                phases[width, metric, comparison, order, q] += 1
        self.assertEqual(len(counts), 1080)
        self.assertEqual(set(counts.values()), {1})
        self.assertEqual(len(phases), 4320)
        self.assertEqual(set(phases.values()), {2})

    def test_control_failure_cannot_be_rescued(self):
        records = copy.deepcopy(self.records)
        for r in records:
            c = r['coordinate']
            if c[5] == 0 and (C.SIDES[c[6]] ^ c[2]) == 1:
                r['observation']['clocks'][3] += 10000
        result = C.statistics(records)
        self.assertEqual(result['outcome'], 'CONTROL_FAIL')
        self.assertFalse(result['WH1_qualified'])
        self.assertFalse(result['all_WH1_faster'])

    def test_wide_decoder_must_improve(self):
        records = copy.deepcopy(self.records)
        for r in records:
            c = r['coordinate']
            if c[3] == 2 and c[4] > 0 and c[5] == 3:
                r['observation']['clocks'][3] = r['observation']['clocks'][2]+100000
        self.assertEqual(C.statistics(records)['outcome'], 'FAIL')

    def test_local_progress_is_not_all_wh1_win(self):
        records = copy.deepcopy(self.records)
        for r in records:
            c = r['coordinate']
            if c[5] == 4 and c[7] == 1:
                r['observation']['clocks'][3] += 50000
        result = C.statistics(records)
        self.assertEqual(result['outcome'], 'PASS')
        self.assertFalse(result['all_WH1_faster'])
        self.assertTrue(result['WH1_qualified'])

    def test_no_sample_omission(self):
        with self.assertRaises(ValueError):
            C.statistics(self.records[:-1])

    def test_clock_corruption_rejected(self):
        for field in ('before', 'after', 'clocks'):
            observation = copy.deepcopy(self.header['prelude'])
            observation[field][0] = True
            with self.assertRaises(ValueError):
                C.clocks(observation, None)

    def test_counter_and_clock_reversal_rejected(self):
        a = copy.deepcopy(self.header['prelude']); a['after'][0] = 2
        with self.assertRaises(ValueError):
            C.clocks(self.header['prelude'], a)
        a = copy.deepcopy(self.header['prelude']); a['clocks'][3] = a['clocks'][2]
        with self.assertRaises(ValueError):
            C.clocks(a, None)

    def test_partial_spool_cleanup(self):
        with tempfile.TemporaryDirectory(prefix='wh2-cost-spool-test-') as d:
            paths = [Path(d)/'first', Path(d)/'second']
            paths[1].write_bytes(b'preserved')
            raw, err, code, failure = C.capture('/not-run', '0'*64, time.monotonic()+1, paths)
            self.assertIsNotNone(failure)
            self.assertEqual((raw, err, code), (b'', b'', None))
            self.assertEqual(paths[0].stat().st_mode & 0o777, 0o400)
            self.assertEqual(paths[1].read_bytes(), b'preserved')

    def test_spent_namespace_no_launch(self):
        with tempfile.TemporaryDirectory(prefix='wh2-cost-spent-test-') as d:
            path = Path(d)/'receipt'
            path.write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL)))
            with patch.object(C, 'OUTPUT', Path(d)), patch.object(C, 'current'), patch.object(C, 'capture') as capture:
                with self.assertRaises(FileExistsError):
                    C.run(path)
                capture.assert_not_called()

    def test_ledger_and_delayed_interference(self):
        records = copy.deepcopy(self.records)
        # A later first callback with a fault remains a valid observation. Shift
        # the complete later chronology; do not subtract it from WORK duration.
        for i, r in enumerate(records):
            shift = (i+1)*1000
            r['ready'] += shift; r['target'] += shift
            r['wait'][0] += shift; r['wait'][2] += shift
            for i in (0, 2, 3, 5):
                r['observation']['clocks'][i] += shift+200
            r['observation']['before'][0] = 0 if r['coordinate'][0] == 0 else 2
            r['observation']['after'][0] = 2
        raw = b''.join(C.A.canonical(r) for r in [self.header]+records+[self.footer])
        self.assertEqual(C.verify(raw, '0'*64, self.old)['outcome'], 'PASS')
        records[-1]['counts'][-1] -= 1
        raw = b''.join(C.A.canonical(r) for r in [self.header]+records+[self.footer])
        with self.assertRaisesRegex(ValueError, 'lifecycle ledger'):
            C.verify(raw, '0'*64, self.old)

    def test_timeout_retains_prefix_and_reaps(self):
        real = subprocess.Popen
        helper = ['/usr/bin/python3', '-c', 'import sys,time; print("prefix",flush=True); time.sleep(2)']
        with tempfile.TemporaryDirectory(prefix='wh2-cost-timeout-test-') as d:
            paths = [Path(d)/'raw', Path(d)/'error']
            with patch.object(C.subprocess, 'Popen', side_effect=lambda args, **kw: real(helper, **kw)):
                raw, error, code, failure = C.capture('/not-a-codec', '0'*64, time.monotonic()+.3, paths)
            self.assertEqual(raw, b'prefix\n')
            self.assertEqual(error, b'')
            self.assertIsNotNone(failure)
            self.assertEqual(code, -9)
            self.assertEqual(paths[0].read_bytes(), raw)

    def test_output_cap_retains_exact_prefix(self):
        real = subprocess.Popen
        helper = ['/usr/bin/python3', '-c', 'print("x"*100,flush=True)']
        with tempfile.TemporaryDirectory(prefix='wh2-cost-cap-test-') as d:
            paths = [Path(d)/'raw', Path(d)/'error']
            with patch.object(C, 'RAW_CAP', 32), patch.object(C.subprocess, 'Popen', side_effect=lambda args, **kw: real(helper, **kw)):
                raw, error, code, failure = C.capture('/not-a-codec', '0'*64, time.monotonic()+2, paths)
            self.assertEqual(raw, b'x'*32)
            self.assertIsNotNone(failure)
            self.assertIsNotNone(code)
            self.assertEqual(paths[0].read_bytes(), raw)


if __name__ == '__main__':
    unittest.main()
