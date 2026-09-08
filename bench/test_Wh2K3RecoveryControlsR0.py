"""Synthetic recovery-receipt checks; never launches a codec workload."""
import copy
import json
import os
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

import Wh2K3RecoveryControlsR0 as C


def synthetic():
    records = []
    for e in C.roster():
        width, ids = e['B'], e['ids']
        arms = []
        for arm in range(3):
            first = C.candidate_first(tuple(ids)) if arm == 1 else (
                0 if e['index'] % (17+arm) == 0 else min(len(ids), 4+arm))
            profile = (struct.pack('<4sHHQQII', b'WHV2', 1, 32, 1, 3*width, width, 0) if arm == 0 else
                       struct.pack('<4sHHQQII', b'WHK3', 1, 32, 0x5748324b33544d31, 3*width, width, 0) if arm == 1 else
                       b'\0'*32)
            packets = [C.candidate_packet(width, i) for i in ids] if arm == 1 else ['00'*32]*len(ids)
            arms.append(dict(profile=profile.hex(), packets=packets,
                             feed=[1]*(first-1)+[0] if first else [1]*len(ids),
                             first=first, recoveries=2 if first else 0, checked=True))
        records.append(dict(type='record', arms=arms, **{k: e[k] for k in ('group', 'index', 'B', 'ids')}))
    return records


def raw(records, mode='native', claim='0'*64):
    header = dict(type='header', protocol=C.PROTOCOL, claim=claim, backend=mode,
                  retained_raw_sha256=C.D.RAW_SHA, features=[0]*4 if mode == 'scalar' else [1]*4,
                  sources=[C.A.sha(C.message(b)) for b in C.WIDTHS])
    footer = dict(type='footer', records=6269, checked=True)
    return b''.join(C.A.canonical(r) for r in [header]+records+[footer])


class RecoveryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = synthetic()
        cls.raw = raw(cls.records)

    def reject(self, mutate, pattern=None):
        records = copy.deepcopy(self.records)
        mutate(records)
        with self.assertRaisesRegex(ValueError, pattern or '.'):
            C.verify(raw(records), '0'*64, 'native')

    def test_complete_synthetic_three_backends(self):
        expected, records = C.verify(self.raw, '0'*64, 'native')
        self.assertEqual(expected['outcome'], 'PASS')
        self.assertTrue(expected['comparative_oh0_win'])
        self.assertFalse(expected['fresh_sample'])
        self.assertFalse(expected['speed_claimed'])
        self.assertEqual([g['cases'] for g in expected['groups']], [6144, 72, 53])
        self.assertEqual(expected['groups'][0]['failure_counts_oh0_to_4'][1], [0]*5)
        self.assertEqual(len(expected['cells']), 12)
        for mode in ('scalar', 'asan'):
            summary, compared = C.verify(raw(self.records, mode), '0'*64, mode)
            self.assertEqual(compared, records)
            self.assertEqual(summary, expected)

    def test_paired_fixes_and_introductions(self):
        result = C.summarize(self.records, C.roster())
        fresh = result['groups'][0]
        for arm, paired in ((0, fresh['paired'][0]), (2, fresh['paired'][1])):
            for overhead in range(5):
                self.assertEqual(len(paired['fixed'][overhead]), fresh['failure_counts_oh0_to_4'][arm][overhead])
                self.assertEqual(paired['introduced'][overhead], [])

    def test_missing_last_case(self):
        with self.assertRaisesRegex(ValueError, 'whole recovery cohort'):
            C.verify(raw(self.records[:-1]), '0'*64, 'native')

    def test_reordered_cases(self):
        def change(rows):
            rows[0], rows[1] = rows[1], rows[0]
        self.reject(change, 'retained exact chronology')

    def test_wrong_width_id_group_or_historical_index(self):
        for key, value in (('B', 64), ('ids', [0,1,2]), ('group', 1), ('index', 5)):
            self.reject(lambda rows: rows[0].__setitem__(key, value), 'retained exact chronology')
        self.reject(lambda rows: rows[-1].__setitem__('index', 0), 'retained exact chronology')

    def test_candidate_packet_corruption(self):
        self.reject(lambda rows: rows[0]['arms'][1]['packets'].__setitem__(0, '00'*32),
                    'independent K3 payload hashes')

    def test_candidate_cannot_inherit_control_failure(self):
        def change(rows):
            arm = rows[0]['arms'][1]
            arm.update(first=0, feed=[1]*7, recoveries=0)
        self.reject(change, 'independent first-success rank')

    def test_every_success_has_two_checked_recoveries(self):
        self.reject(lambda rows: rows[1]['arms'][2].__setitem__('recoveries', 1), 'two byte-exact recoveries')
        self.reject(lambda rows: rows[-1]['arms'][0].__setitem__('checked', False), 'payload/guard validation')

    def test_control_unresolved_is_retained_not_rejected(self):
        self.assertEqual(self.records[0]['arms'][0]['first'], 0)
        summary, _ = C.verify(self.raw, '0'*64, 'native')
        self.assertGreater(summary['groups'][0]['unresolved'][0], 0)

    def test_feed_endpoint_and_prefix_must_match(self):
        for first in (1, 2, 8, True):
            self.reject(lambda rows: rows[0]['arms'][0].__setitem__('first', first))
        self.reject(lambda rows: rows[0]['arms'][0]['feed'].__setitem__(0, 0), 'all prefix feed statuses')
        self.reject(lambda rows: rows[1]['arms'][0]['feed'].append(0), 'all prefix feed statuses')

    def test_descriptor_identity_and_stability(self):
        for arm in (0, 1, 2):
            self.reject(lambda rows: rows[0]['arms'][arm].__setitem__('profile', 'ff'*32))
        self.reject(lambda rows: rows[-1]['arms'][0].__setitem__('profile', '00'*32),
                    'stable full descriptor')

    def test_bad_hash_type_or_packet_omission(self):
        for value in ('00', 'GG'*32, True, None):
            self.reject(lambda rows: rows[0]['arms'][0]['packets'].__setitem__(0, value), 'hex bytes')
        self.reject(lambda rows: rows[-1]['arms'][2]['packets'].pop(), 'all retained packets')

    def test_truncation_and_claim_backend_sealed_source(self):
        for broken in (self.raw[:-1], self.raw[:100], self.raw.replace(b'"checked":true', b'"checked":false', 1)):
            with self.assertRaises(ValueError):
                C.verify(broken, '0'*64, 'native')
        for claim, mode in (('1'*64, 'native'), ('0'*64, 'scalar')):
            with self.assertRaisesRegex(ValueError, 'exact replay header'):
                C.verify(self.raw, claim, mode)
        broken = self.raw.replace(C.D.RAW_SHA.encode(), b'0'*64, 1)
        with self.assertRaisesRegex(ValueError, 'exact replay header'):
            C.verify(broken, '0'*64, 'native')

    def test_spent_namespace_never_launches(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k3-recovery-spent-') as d:
            path = Path(d)/'receipt'
            path.write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL)))
            with patch.object(C, 'OUTPUT', Path(d)), patch.object(C, 'current'), patch.object(C.U, 'capture') as capture:
                with self.assertRaises(FileExistsError):
                    C.run(path)
                capture.assert_not_called()

    def test_receipt_cannot_drop_pins_or_rebind_executable(self):
        files, manifests, executables = {}, {}, {}
        def add(path, content):
            files[path] = dict(path=path, bytes=len(content), sha256=C.A.sha(content))
            return files[path]
        source = add('/synthetic/source.cpp', b'source')
        for mode in C.MODES:
            executable = '/synthetic/'+mode+'/recovery_worker'
            executables[mode] = executable
            artifact = add(executable, mode.encode())
            path = str(Path(executable).parent/'manifest.json')
            manifests[path] = C.A.canonical(dict(protocol=C.PROTOCOL, mode=mode,
                                                inputs=[source], artifacts=[artifact]))
            add(path, manifests[path])
        environment = dict({k: None for k in C.U.ENV_KEYS}, **C.SANITIZERS)
        receipt = dict(protocol=C.PROTOCOL, head='head', executables=executables,
                       environment=environment, pins=sorted(files.values(), key=lambda p: p['path']))
        def read(path, cap):
            return manifests[str(path)]
        with patch.dict(os.environ, C.SANITIZERS), patch.object(C.U, 'command', return_value=b'head\n'), \
                patch.object(C.U, 'pin', side_effect=lambda p: files[str(p)]), \
                patch.object(C.A, 'read_regular', side_effect=read), patch.object(C, 'retained'):
            C.current(receipt)
            for mutate in (lambda r: r['pins'].remove(source),
                           lambda r: r['pins'].append(source),
                           lambda r: r['executables'].__setitem__('native', '/unbound/native/recovery_worker'),
                           lambda r: r['executables'].__setitem__('native', executables['scalar'])):
                changed = copy.deepcopy(receipt); mutate(changed)
                with self.assertRaises(ValueError):
                    C.current(changed)

    def test_cross_backend_packet_mismatch_is_invalid(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k3-recovery-parity-') as d:
            path = Path(d)/'receipt'
            path.write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL, executables={m: m for m in C.MODES})))
            def capture(executable, claim, deadline, spools):
                records = copy.deepcopy(self.records)
                if executable == 'scalar':
                    records[0]['arms'][0]['packets'][0] = '11'*32
                data = raw(records, executable, claim)
                C.A.publish(spools[0], data); C.A.publish(spools[1], b'')
                return data, b'', 0, None
            output = Path(d)/'result'
            with patch.object(C, 'OUTPUT', output), patch.object(C, 'current'), patch.object(C.U, 'capture', side_effect=capture):
                C.run(path)
            analysis = json.loads((output/'analysis.json').read_bytes())
            self.assertEqual(analysis['outcome'], 'INVALID')
            self.assertIn('cross-backend complete byte/status parity', analysis['failure'])
            self.assertFalse((output/'asan.raw.jsonl').exists())
            self.assertTrue((output/'COMPLETE.json').exists())


if __name__ == '__main__':
    unittest.main()
