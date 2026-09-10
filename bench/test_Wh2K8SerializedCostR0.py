"""Synthetic timing/ledger tests only: never invokes a scientific codec worker."""
import copy
import json
import os
from pathlib import Path
import tempfile
import time
import unittest
import subprocess
from unittest.mock import patch

import Wh2K8SerializedCostR0 as C
import struct


def synthetic():
    fixtures, old = [], []
    rows = bytes(v for row in C.selected_rows() for v in row).hex()
    for b,tail in C.SHAPES:
        candidate_profile, candidate_packets = C.candidate_fixture(b,tail)
        generic = dict(profile='00'*32, packets=candidate_packets, steps=[C.K,C.K], rows=rows, fallback=[False,False])
        candidate = dict(profile=candidate_profile, packets=candidate_packets, steps=[C.K,C.K], rows=rows, fallback=[False,False])
        certified = copy.deepcopy(generic)
        certified['profile'] = struct.pack('<4sHHQQII',b'WHV2',1,32,0x4b295bbb47f4f9c9,(C.K-1)*b+tail,b,0).hex()
        arms = [certified,candidate,generic]
        fixtures.append(dict(width=b,message=(C.K-1)*b+tail,packet_bytes=C.packet_lengths(b,tail),source=bytes((37*i+i//11) % 256 for i in range((C.K-1)*b+tail)).hex(),
                             arms=arms+copy.deepcopy(arms)))
        old.append(copy.deepcopy(fixtures[-1]))
    prelude = dict(clocks=[100, 100, 102, 112, 108, 120], before=[0]*4, after=[0]*4)
    header = dict(type='header', protocol=C.PROTOCOL, claim='0'*64, batch=C.BATCH,
                  identity_hex='00', prelude=prelude, fixtures=fixtures)
    records = []; previous = prelude; total = 0
    for coordinate in C.roster():
        metric, arm, q = coordinate[4], coordinate[7], coordinate[8]
        ready, cpu = previous['clocks'][5]+100, previous['clocks'][4]
        target = ready+q; duration = (100000,98000,120000,100000,97000,120000)[arm]
        observation = dict(clocks=[target+4, cpu+40, target+6, target+6+duration,
                                  cpu+40+duration, target+20+duration], before=[0]*4, after=[0]*4)
        counts = [0 if metric else 128, 0 if metric else C.BATCH*len(C.PACKETS), 128 if metric else 0,
                  C.BATCH*C.K if metric else 0, 128 if metric else 0, 128]
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
        self.assertTrue(result['speed_qualified'])
        self.assertEqual(len(result['statistics']), 360)
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
        self.assertEqual(len(counts), 4320)
        self.assertEqual(set(counts.values()), {1})
        self.assertEqual(len(phases), 17280)
        self.assertEqual(set(phases.values()), {2})

    def test_control_failure_cannot_be_rescued(self):
        records = copy.deepcopy(self.records)
        for r in records:
            c = r['coordinate']
            if c[5] == 0 and (C.SIDES[c[6]] ^ c[2]) == 1:
                r['observation']['clocks'][3] += 10000
        result = C.statistics(records)
        self.assertEqual(result['outcome'], 'CONTROL_FAIL')
        self.assertFalse(result['speed_qualified'])

    def test_every_control_cell_enforces_both_strict_bounds(self):
        bound=C.math.log1p(.02); checked=0
        for width in range(6):
            for metric in range(3):
                for comparison in range(6):
                    for order in range(2):
                        key=(width,metric,comparison,order)
                        for low,high,passed in ((-bound,0,False),(0,bound,False),
                                (-bound-1e-12,0,False),(0,bound+1e-12,False),
                                (-bound+1e-12,bound-1e-12,True),(0,0,True)):
                            with patch.object(C.A,'confidence',return_value=dict(lower95_log=low,upper95_log=high)):
                                item=C.score_cell(key,[0]*12)
                            self.assertEqual(item['control_pass'],passed)
                            self.assertNotIn('treatment_pass',item)
                            self.assertEqual(item['comparison_arms'],[C.ARM_NAMES[comparison]]*2)
                        checked+=1
        self.assertEqual(checked,216)

    def test_wide_decoder_must_improve(self):
        records = copy.deepcopy(self.records)
        for r in records:
            c = r['coordinate']
            if c[3] == 4 and c[4] > 0 and c[5] == 6:
                r['observation']['clocks'][3] = r['observation']['clocks'][2]+100000
        self.assertEqual(C.statistics(records)['outcome'], 'FAIL')

    def test_must_beat_wh1_too(self):
        records = copy.deepcopy(self.records)
        for r in records:
            c = r['coordinate']
            if c[5] == 7 and c[7] == 1:
                r['observation']['clocks'][3] += 50000
        result = C.statistics(records)
        self.assertEqual(result['outcome'], 'FAIL')
        self.assertFalse(result['speed_qualified'])

    def test_every_treatment_cell_requires_improvement(self):
        # A tiny regression/equality in ANY width, lifecycle, order or control
        # is not a pass; there is no legacy two-percent treatment tolerance.
        for width in range(6):
            for metric in range(3):
                for comparison in (6, 7, 8, 9):
                    for order in range(2):
                        records = [dict(coordinate=r['coordinate'], observation=dict(clocks=list(r['observation']['clocks'])))
                                   for r in self.records]
                        for r in records:
                            c = r['coordinate']
                            if (c[3], c[4], c[5], c[2]) == (width, metric, comparison, order):
                                r['observation']['clocks'][3] = r['observation']['clocks'][2]+100000
                        result = C.statistics(records)
                        self.assertEqual(result['outcome'], 'FAIL')
                        self.assertIn([width,metric,comparison,order], result['failed_treatments'])

    def test_candidate_oracle_and_descriptor_are_required(self):
        header = copy.deepcopy(self.header)
        header['fixtures'][0]['arms'][1]['profile'] = '00'*32
        with self.assertRaises(ValueError):
            C.verify_fixtures(header, self.old)
        header = copy.deepcopy(self.header)
        header['fixtures'][0]['arms'][1]['packets'] = '00'*(2*len(C.PACKETS))
        header['fixtures'][0]['arms'][4]['packets'] = '00'*(2*len(C.PACKETS))
        old = copy.deepcopy(self.old)
        old['fixtures'][0]['arms'][1]['packets'] = '00'*(2*len(C.PACKETS))
        old['fixtures'][0]['arms'][4]['packets'] = '00'*(2*len(C.PACKETS))
        with self.assertRaisesRegex(ValueError, 'independent native-row packet arithmetic'):
            C.verify_fixtures(header, old)

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

    def test_invalid_partial_spool_bundle_replays(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k8-cost-invalid-test-') as d:
            base=Path(d); output=base/'scientific'; receipt=base/'receipt'
            receipt.write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL)))
            def fail(executable,claim,deadline,spools):
                del executable,claim,deadline
                spools[0].write_bytes(b'partial\n')
                spools[0].chmod(0o400)
                return b'partial\n',b'',None,'second spool creation failed'
            receipt.write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL,executable='/not-run')))
            with patch.object(C,'OUTPUT',output),patch.object(C,'current'),patch.object(C,'capture',side_effect=fail):
                C.run(receipt)
                self.assertEqual(C.replay()['outcome'],'INVALID')
                self.assertFalse((output/'stderr.txt').exists())

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


    def test_binary_target_roundtrip(self):
        binary = b'!' + b'\x00' + bytes(range(256))*2 + bytes(range(103))
        self.assertEqual(len(binary), 617)
        target = copy.deepcopy(self.header)
        target['identity_hex'] = binary.hex()
        def read(path, cap):
            return C.A.canonical(target if path.name == 'target.json' else self.header)
        with patch.object(C.A, 'read_regular', side_effect=read), patch.object(C, 'TARGET_SHA', C.A.sha(binary)):
            result = C.prior_header(Path('/synthetic-not-read'))
        self.assertEqual(bytes.fromhex(result['identity_before']['canonical_hex']), binary)

    def test_truncated_corrupt_target_rejected_before_launch(self):
        binary = b'!' + b'\x00' + bytes(range(256))*2 + bytes(range(103))
        for encoded in ('21', binary[:-1].hex(), binary.hex()[:-1], binary.hex()+'00',
                        binary.hex().upper(), 'g'*1234, '00'*617, None):
            target = copy.deepcopy(self.header)
            target['identity_hex'] = encoded
            def read(path, cap):
                return C.A.canonical(target if path.name == 'target.json' else self.header)
            with patch.object(C.A, 'read_regular', side_effect=read), patch.object(C, 'TARGET_SHA', C.A.sha(binary)):
                with self.assertRaises(ValueError):
                    C.prior_header(Path('/synthetic-not-read'))

    def test_target_fixture_and_schema_must_match(self):
        binary = b'!' + b'\x00' + bytes(range(256))*2 + bytes(range(103))
        for change in (lambda h: h.update(protocol='wrong'), lambda h: h.update(batch=True),
                       lambda h: h['fixtures'][0]['arms'][0]['steps'].__setitem__(0,6),
                       lambda h: h.pop('prelude')):
            target = copy.deepcopy(self.header); target['identity_hex'] = binary.hex()
            change(target)
            def read(path, cap):
                return C.A.canonical(target if path.name == 'target.json' else self.header)
            with patch.object(C.A, 'read_regular', side_effect=read), patch.object(C, 'TARGET_SHA', C.A.sha(binary)):
                with self.assertRaises(ValueError):
                    C.prior_header(Path('/synthetic-not-read'))

    def test_original_percent_s_output_is_not_a_reference(self):
        def read(path, cap):
            return b'!\n' if path.name == 'target.json' else C.A.canonical(self.header)
        with patch.object(C.A, 'read_regular', side_effect=read):
            with self.assertRaises(ValueError):
                C.prior_header(Path('/synthetic-not-read'))

    def test_control_cannot_silently_use_the_new_default(self):
        header = copy.deepcopy(self.header)
        header['fixtures'][0]['arms'][0]['profile'] = header['fixtures'][0]['arms'][1]['profile']
        old = copy.deepcopy(self.old)
        old['fixtures'] = copy.deepcopy(header['fixtures'])
        with self.assertRaisesRegex(ValueError, 'explicit certified WH2 descriptor'):
            C.verify_fixtures(header, old)

    def test_both_candidate_policies_use_identical_equations(self):
        for field in ('profile', 'packets', 'steps', 'rows'):
            header = copy.deepcopy(self.header)
            header['fixtures'][0]['arms'][4][field] = header['fixtures'][0]['arms'][0][field]
            if field in ('packets','rows'):
                header['fixtures'][0]['arms'][4][field] = 'ff'+header['fixtures'][0]['arms'][4][field][2:]
            if field == 'steps':
                header['fixtures'][0]['arms'][4][field] = [6, 5]
            old = copy.deepcopy(self.old)
            old['fixtures'] = copy.deepcopy(header['fixtures'])
            with self.assertRaises(ValueError):
                C.verify_fixtures(header, old)

    def test_receipt_cannot_drop_pins_or_rebind_executable(self):
        files = {}
        def add(path, content):
            files[path] = dict(path=path, bytes=len(content), sha256=C.A.sha(content))
            return files[path]
        source = add('/synthetic/source.cpp', b'source')
        executable = '/synthetic/native/cost_worker'
        artifact = add(executable, b'worker')
        manifest_path = '/synthetic/native/manifest.json'
        manifest = C.A.canonical(dict(protocol=C.PROTOCOL, mode='native',
                                     inputs=[source], artifacts=[artifact]))
        add(manifest_path, manifest)
        receipt = dict(protocol=C.PROTOCOL, head='head', executable=executable,
                       environment={k: None for k in C.ENV_KEYS},
                       pins=sorted(files.values(), key=lambda p: p['path']))
        with patch.dict(os.environ, {}, clear=True), patch.object(C, 'command', return_value=b'head\n'), \
                patch.object(C, 'pin', side_effect=lambda p: files[str(p)]), \
                patch.object(C.A, 'read_regular', return_value=manifest), patch.object(C, 'prior_header'):
            C.current(receipt)
            for mutate in (lambda r: r['pins'].remove(source), lambda r: r['pins'].append(source),
                           lambda r: r.__setitem__('executable', '/unbound/native/cost_worker'),
                           lambda r: r.__setitem__('executable', '/synthetic/scalar/cost_worker'),
                           lambda r: r.__setitem__('executable', source['path'])):
                changed = copy.deepcopy(receipt); mutate(changed)
                with self.assertRaises(ValueError):
                    C.current(changed)


    def test_control_rows_and_payloads_are_not_trusted(self):
        for a in (0,2,3,5):
            for key in ('packets','rows'):
                header = copy.deepcopy(self.header)
                header['fixtures'][0]['arms'][a][key] = 'ff'+header['fixtures'][0]['arms'][a][key][2:]
                prior = copy.deepcopy(self.old); prior['fixtures'] = copy.deepcopy(header['fixtures'])
                with self.subTest(arm=a,key=key), self.assertRaises(ValueError):
                    C.verify_fixtures(header,prior)

    def test_required_frozen_sizes_and_ownership_pairs(self):
        self.assertEqual(C.K,8)
        self.assertEqual(C.PAIRS,((0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,1),(2,1),(3,4),(5,4)))
        self.assertEqual(C.CALLBACKS,77760)
        self.assertEqual(C.ARM_NAMES,('ordinary_independent','k8_independent','wh1_owned',
                                     'ordinary_borrowed','k8_borrowed','wh1_borrowed'))
        self.assertEqual(len(C.selected_rows()),24)
        self.assertEqual(C.RAW_CAP,192*1024**2)
        self.assertNotEqual(C.OUTPUT,Path('/var/tmp/wh2-k3-ordinary-cost-r0'))

    def publication(self,mode):
        header=copy.deepcopy(self.header); header['identity_hex']=b'neutral-deferred'.hex()
        records=copy.deepcopy([self.records[i] for i in (0,180,360)])
        for row in records:
            row.update(ready=0,target=0,wait=[0]*4)
        last=records[-1]
        if mode!='success': last['checked']=False
        if mode in ('last-recover','throw-recover'): last['complete']=False
        if mode=='last-clock': last['observation']['clocks'][5]=0
        total=sum(r['observation']['clocks'][3]-r['observation']['clocks'][2]
                  for r in (records[:-1] if mode=='last-clock' else records))
        footer=dict(type='footer',complete=mode=='success',records=3,work_ns=total)
        return [header]+records+[footer]

    def test_neutral_publication_and_late_failures(self):
        for mode in ('success','last-recover','throw-recover','last-clock','last-source'):
            rows=self.publication(mode)
            raw=b''.join(C.A.canonical(r) for r in rows)
            self.assertEqual(C.verify_publication(raw,mode,self.old),
                             dict(neutral=True,records=3,mode=mode,speed_qualified=False))
            with self.assertRaises(ValueError): C.verify(raw,'0'*64,self.old)

    def test_publication_rejects_omission_retry_and_corruption(self):
        for mode in ('success','last-recover','throw-recover','last-clock','last-source'):
            for mutation in (
                    lambda r:r.pop(3), lambda r:r.insert(1,copy.deepcopy(r[0])),
                    lambda r:r[3]['counts'].__setitem__(5,127),
                    lambda r:r[3]['addresses'].__setitem__(127,0),
                    lambda r:r[3].__setitem__('checked',not r[3]['checked']),
                    lambda r:r[3].__setitem__('complete',not r[3]['complete']),
                    lambda r:r[-1].__setitem__('work_ns',r[-1]['work_ns']+1),
                    lambda r:r[-1].__setitem__('complete',not r[-1]['complete']),
                    lambda r:r[3]['coordinate'].__setitem__(0,359)):
                rows=self.publication(mode); mutation(rows)
                with self.subTest(mode=mode),self.assertRaises(ValueError):
                    C.verify_publication(b''.join(C.A.canonical(r) for r in rows),mode,self.old)

    def test_old_scientific_namespace_and_work_unchanged(self):
        old=(C.ROOT/'bench/Wh2K5DeferredCostR0.cpp').read_text()
        new=(C.ROOT/'bench/Wh2SmallLifecycleWorkerR0.h').read_text()
        for start,end in (('struct Reader {','struct Work {'),):
            self.assertEqual(old[old.index(start):old.index(end)],new[new.index(start):new.index(end)])
        loop=new[new.index('int Worker('):new.index('struct FakeReader')]
        measured=loop[:loop.index('} catch(const std::exception&')]
        self.assertNotIn('RecordJson(',measured); self.assertNotIn('HeaderJson(',measured)
        self.assertNotIn('Publish(',measured)
        self.assertEqual(loop.count('Publish('),1)
        self.assertEqual(C.OUTPUT,Path('/var/tmp/wh2-k8-serialized-cost-r0'))

    def test_partial_tail_bytes_and_guards(self):
        for shape,(b,tail) in enumerate(C.SHAPES):
            f=self.header['fixtures'][shape]
            self.assertEqual(f['message'],(C.K-1)*b+tail)
            self.assertEqual(f['packet_bytes'],[b]*7+[tail]+[b]*16)
            packets=bytes.fromhex(f['arms'][1]['packets'])
            self.assertEqual(packets[7*b+tail:8*b],bytes([0xa5])*(b-tail))
            if tail<b:
                changed=copy.deepcopy(self.header)
                changed['fixtures'][shape]['arms'][1]['packets']=(
                    packets[:7*b+tail]+b'\x00'+packets[7*b+tail+1:]).hex()
                prior=copy.deepcopy(self.old);prior['fixtures']=copy.deepcopy(changed['fixtures'])
                with self.assertRaises(ValueError):C.verify_fixtures(changed,prior)

    def test_actual_feed_endpoint_and_fallback_are_charged(self):
        header=copy.deepcopy(self.header);f=header['fixtures'][1];b,tail=C.SHAPES[1]
        # Zero control repair rows require all eight systematic fallback packets.
        # Never equalize control work to the candidate endpoint.
        rows=tuple(tuple(int(i==j) for j in range(C.K)) for i in range(C.K))+((0,)*C.K,)*(2*C.REPAIRS)
        for arm in (0,3):
            f['arms'][arm].update(rows=bytes(v for row in rows for v in row).hex(),
                packets=C.payload(rows,b,tail),steps=[16,16],fallback=[True,True])
        prior=copy.deepcopy(self.old);prior['fixtures']=copy.deepcopy(header['fixtures'])
        C.verify_fixtures(header,prior)
        records=copy.deepcopy(self.records)
        for record in records:
            c=record['coordinate']
            if c[3]==1 and c[4]>0 and c[7] in (0,3):record['counts'][3]=128*16
        raw=b''.join(C.A.canonical(row) for row in [header]+records+[self.footer])
        self.assertEqual(C.verify(raw,'0'*64,prior)['outcome'],'PASS')
        wrong=copy.deepcopy(header);wrong['fixtures'][1]['arms'][0]['fallback']=[False,False]
        wrong_prior=copy.deepcopy(prior);wrong_prior['fixtures']=copy.deepcopy(wrong['fixtures'])
        with self.assertRaises(ValueError):C.verify_fixtures(wrong,wrong_prior)
        records[0]['counts'][1]-=1
        raw=b''.join(C.A.canonical(row) for row in [header]+records+[self.footer])
        with self.assertRaises(ValueError):C.verify(raw,'0'*64,prior)



if __name__ == '__main__':
    unittest.main()
