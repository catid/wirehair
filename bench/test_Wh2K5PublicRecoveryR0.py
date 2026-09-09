"""Synthetic recovery and controller tests; never launches a codec worker."""
import copy
import json
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import time
import unittest
from unittest.mock import patch

import Wh2K5PublicRecoveryR0 as C


def arm_result(width, ids, arm, deficient=False):
    rows = tuple((1,0,0,0,0) if deficient else C.coefficient(i) for i in ids)
    first = C.first_success(rows)
    if arm in (0,3):
        profile = struct.pack('<4sHHQQII',b'WHV2',1,32,0x4b295bbb47f4f9c9,5*width,width,0)
    elif arm in (1,4):
        profile = struct.pack('<4sHHQQII',b'WHV2',1,32,0x80070c81bfe375f1,5*width,width,0)
    else: profile = bytes(32)
    return dict(profile=profile.hex(),rows=bytes(v for row in rows for v in row).hex(),
                packets=[C.packet_hash(width,row) for row in rows],
                feed=[1]*(first-1)+[0] if first else [1]*len(ids),first=first,
                recoveries=2 if first else 0,counts=[6,6*len(ids),1,first or len(ids),2 if first else 0,7],checked=True)


def synthetic(neutral=False):
    records = []
    for e in C.roster(neutral):
        arms = [arm_result(e['B'],e['ids'],a,deficient=a!=1 and e['index'] % (17+a)==0) for a in range(3)]
        records.append(dict(type='record',arms=arms+copy.deepcopy(arms),
                            **{k:e[k] for k in ('group','index','B','ids')}))
    return records


def raw(records, mode='native', claim='0'*64, neutral=False):
    header = dict(type='header',protocol=C.PROTOCOL,claim=claim,backend=mode,
                  scope='neutral' if neutral else 'retained',retained_raw_sha256=C.D.RAW_SHA,
                  features=[0]*4 if mode=='scalar' else [1]*4,
                  sources=[C.A.sha(C.message(b)) for b in C.WIDTHS])
    return b''.join(C.A.canonical(r) for r in [header]+records+
                    [dict(type='footer',records=36 if neutral else 6273,checked=True)])


class RecoveryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = synthetic()
        cls.raw = raw(cls.records)

    def reject(self, change, pattern=None):
        records = copy.deepcopy(self.records); change(records)
        with self.assertRaisesRegex(ValueError,pattern or '.'):
            C.verify(raw(records),'0'*64,'native')

    def test_complete_three_backends(self):
        summary,records = C.verify(self.raw,'0'*64,'native')
        self.assertEqual(summary['outcome'],'PASS')
        self.assertTrue(summary['comparative_oh0_win'])
        self.assertTrue(summary['candidate_cells_at_most_one_percent'])
        self.assertEqual([g['cases'] for g in summary['groups']],[6144,72,57])
        for a in (1,4): self.assertEqual(summary['groups'][0]['failure_counts_oh0_to_4'][a],[11,0,0,0,0])
        for flag in ('fresh_sample','speed_claimed','all_K_claimed','production_promotion_claimed'):
            self.assertFalse(summary[flag])
        self.assertEqual(len(summary['cells']),12)
        for mode in ('scalar','asan'):
            other,compared = C.verify(raw(records,mode),'0'*64,mode)
            self.assertEqual(other,summary); self.assertEqual(compared,records)

    def test_neutral_bounded_smoke(self):
        records = synthetic(True)
        summary,_ = C.verify(raw(records,neutral=True),'0'*64,'native',True)
        self.assertEqual(summary,dict(outcome='PASS',neutral_cases=36))
        smoke = C.roster(True)[12:]
        self.assertEqual([r['index'] for r in smoke],list(range(12))+list(range(6132,6144)))
        cells = {}
        for e in smoke:
            row = C.retained()[0]['fresh'][e['index']]
            key = row['B'],row['schedule']; cells[key] = cells.get(key,0)+1
        self.assertEqual(len(cells),12); self.assertEqual(set(cells.values()),{2})

    def test_coefficients_match_all_sealed_rows(self):
        for r in C.retained()[1]['rows']:
            self.assertEqual(C.coefficient(r['id']),tuple(r['row']))

    def test_paired_differences_and_cells(self):
        summary = C.summarize(self.records,C.roster())
        for g in summary['groups'][:2]+summary['cells']:
            counts = g['failure_counts_oh0_to_4']
            for p,(candidate,control) in zip(g['paired'],C.PAIRED):
                self.assertEqual((p['candidate'],p['control']),(C.ARMS[candidate],C.ARMS[control]))
                for oh in range(5):
                    self.assertEqual(len(p['fixed'][oh])-len(p['introduced'][oh]),counts[control][oh]-counts[candidate][oh])

    def test_introductions_count_in_both_policies(self):
        records = copy.deepcopy(self.records)
        for a in (0,2,3,5): records[0]['arms'][a]['first'] = 5
        for a in (1,4): records[0]['arms'][a]['first'] = 0
        summary = C.summarize(records,C.roster())
        key = [records[0]['index'],records[0]['B']]
        for p in summary['groups'][0]['paired']:
            for oh in range(5):
                self.assertIn(key,p['introduced'][oh]); self.assertNotIn(key,p['fixed'][oh])

    def test_whole_shape_and_chronology(self):
        for change in (lambda r:r.pop(),lambda r:r[0]['arms'].pop(),
                       lambda r:r[0].__setitem__('index',1),lambda r:r[0].__setitem__('B',64),
                       lambda r:r[-1]['ids'].append(0)):
            self.reject(change)

    def test_control_cannot_select_candidate(self):
        self.reject(lambda r:r[0]['arms'][0].__setitem__('profile',r[0]['arms'][1]['profile']),
                    'explicit certified WH2 descriptor')

    def test_every_arm_packet_rows_and_ledger_checked(self):
        for a in range(6):
            self.reject(lambda r:r[0]['arms'][a]['packets'].__setitem__(0,'00'*32),'independent every-arm packet hashes')
            self.reject(lambda r:r[0]['arms'][a]['counts'].__setitem__(0,5),'whole attempted API ledger')
            self.reject(lambda r:r[0]['arms'][a].__setitem__('rows','ff'+r[0]['arms'][a]['rows'][2:]))

    def test_rank_endpoint_and_feed_checked_for_controls(self):
        # Valid source-equation rank forbids a control's fabricated late success.
        self.reject(lambda r:r[1]['arms'][2].__setitem__('first',6),'independent every-arm first-success rank')
        self.reject(lambda r:r[1]['arms'][2]['feed'].__setitem__(0,0),'all prefix feed statuses')

    def test_candidate_preserves_recorded_deficiencies(self):
        index = next(i for i,r in enumerate(self.records[:6144]) if r['arms'][1]['first']==6)
        self.reject(lambda r:r[index]['arms'][4].__setitem__('first',5),'independent every-arm first-success rank')

    def test_unresolved_control_is_retained(self):
        summary,_ = C.verify(self.raw,'0'*64,'native')
        self.assertGreater(summary['groups'][0]['unresolved'][0],0)
        self.assertEqual(summary['groups'][0]['unresolved'][1],0)

    def test_exact_recovery_and_boolean_types(self):
        for a in range(6):
            self.reject(lambda r:r[1]['arms'][a].__setitem__('recoveries',1),'two byte-exact recoveries')
            self.reject(lambda r:r[1]['arms'][a].__setitem__('checked',False),'payload/guard validation')
        for value in (True,1,4,10):
            self.reject(lambda r:r[0]['arms'][0].__setitem__('first',value))

    def test_hex_bounds_and_omitted_packet(self):
        for value in ('00','GG'*32,True,None):
            self.reject(lambda r:r[0]['arms'][0]['packets'].__setitem__(0,value),'hex bytes')
        self.reject(lambda r:r[-1]['arms'][2]['packets'].pop(),'all retained packets')

    def test_claim_backend_truncation_and_scope(self):
        for broken in (self.raw[:-1],self.raw[:100],self.raw.replace(C.D.RAW_SHA.encode(),b'0'*64,1)):
            with self.assertRaises(ValueError): C.verify(broken,'0'*64,'native')
        for claim,mode in (('1'*64,'native'),('0'*64,'scalar')):
            with self.assertRaisesRegex(ValueError,'exact replay header'): C.verify(self.raw,claim,mode)
        with self.assertRaises(ValueError): C.verify(self.raw,'0'*64,'native',True)

    def test_spent_namespace_never_launches(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k5-recovery-spent-') as d:
            receipt = Path(d)/'receipt'; receipt.write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL)))
            with patch.object(C,'OUTPUT',Path(d)),patch.object(C,'current'),patch.object(C,'capture') as capture:
                with self.assertRaises(FileExistsError): C.run(receipt)
                capture.assert_not_called()

    def test_prototype_descriptor_and_seed_rejected(self):
        for profile in (struct.pack('<4sHHQQII',b'WHK5',1,32,0x5748324b35544d31,10,2,0),
                        struct.pack('<4sHHQQII',b'WHV2',1,32,0x80070c81bfe375f1,10,2,1)):
            self.reject(lambda r:r[0]['arms'][1].__setitem__('profile',profile.hex()),
                        'installed sealed WHV2 K5 identity')

    def test_context_symbol_must_be_unique_and_sized(self):
        for size in (0x22810,0x21810,0x26820):
            self.assertEqual(C.context_size(('00000000 %08x B GF256Ctx\n' % size).encode()),size)
        for broken in (b'',b'00000000 B GF256Ctx\n',b'0 nope B GF256Ctx\n',
                       b'0 00000000 B GF256Ctx\n',b'0 00100000 B GF256Ctx\n',
                       b'0 22810 T GF256Ctx\n',b'0 22810 B GF256Ctx\n0 22810 B GF256Ctx\n'):
            with self.assertRaises(ValueError): C.context_size(broken)

    def test_whole_prototype_parity_only_rebinds_candidate_descriptor(self):
        records = copy.deepcopy(self.records)
        for row in records:
            for arm in (1,4):
                row['arms'][arm]['profile'] = struct.pack('<4sHHQQII',b'WHK5',1,32,
                    0x5748324b35544d31,5*row['B'],row['B'],0).hex()
        data = raw(records)
        with patch.object(C,'PROTOTYPE_SHA',C.A.sha(data)),patch.object(C.A,'read_regular',return_value=data):
            C.retained_parity(self.records)
            changed = copy.deepcopy(self.records)
            changed[-1]['arms'][2]['first'] = 0
            with self.assertRaisesRegex(ValueError,'complete installed versus retained prototype parity'):
                C.retained_parity(changed)
        with patch.object(C.A,'read_regular',return_value=b'wrong'):
            with self.assertRaisesRegex(ValueError,'authenticated prototype recovery evidence'):
                C.retained_parity(self.records)

    def test_partial_spool_cleanup(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k5-recovery-spool-') as d:
            paths = [Path(d)/'raw',Path(d)/'error']; paths[1].write_bytes(b'preserved')
            raw_bytes,error,code,failure = C.capture('/never-run','0'*64,time.monotonic()+1,paths)
            self.assertEqual((raw_bytes,error,code),(b'',b'',None)); self.assertIsNotNone(failure)
            self.assertEqual(paths[0].stat().st_mode & 0o777,0o400)
            self.assertEqual(paths[1].read_bytes(),b'preserved')

    def test_receipt_cannot_drop_pins_or_rebind_a_backend(self):
        files,manifests,executables = {},{},{}
        def add(path,content):
            p = dict(path=path,bytes=len(content),sha256=C.A.sha(content))
            files[path] = p; return p
        source = add('/synthetic/source.cpp',b'source')
        for mode in C.MODES:
            executable = '/synthetic/'+mode+'/recovery_worker'; executables[mode] = executable
            artifact = add(executable,mode.encode())
            fixture = add('/synthetic/'+mode+'/fixtures.jsonl',b'fixture')
            path = '/synthetic/'+mode+'/manifest.json'
            manifests[path] = C.A.canonical(dict(protocol=C.PROTOCOL,mode=mode,
                                                inputs=[source],artifacts=[artifact,fixture]))
            add(path,manifests[path])
        environment = dict({k:None for k in C.U.ENV_KEYS},**C.SANITIZERS)
        receipt = dict(protocol=C.PROTOCOL,head='head',executables=executables,
                       environment=environment,pins=sorted(files.values(),key=lambda p:p['path']))
        def read(path,cap):
            return b'fixture' if Path(path).name=='fixtures.jsonl' else manifests[str(path)]
        with patch.dict(os.environ,C.SANITIZERS,clear=True),patch.object(C.U,'command',return_value=b'head\n'), \
                patch.object(C.U,'pin',side_effect=lambda p:files[str(p)]), \
                patch.object(C.A,'read_regular',side_effect=read),patch.object(C,'retained'), \
                patch.object(C,'verify',return_value=({},[])):
            C.current(receipt)
            for mutate in (lambda r:r['pins'].remove(source),lambda r:r['pins'].append(source),
                           lambda r:r['executables'].__setitem__('native','/unbound/native/recovery_worker'),
                           lambda r:r['executables'].__setitem__('native',executables['scalar']),
                           lambda r:r['environment'].__setitem__('ASAN_OPTIONS','detect_leaks=0')):
                changed = copy.deepcopy(receipt); mutate(changed)
                with self.assertRaises(ValueError): C.current(changed)

    def test_output_cap_and_timeout_reap(self):
        real = subprocess.Popen
        for helper,cap,deadline in ((['/usr/bin/python3','-c','print("x"*100,flush=True)'],32,2),
                (['/usr/bin/python3','-c','import time; print("prefix",flush=True); time.sleep(2)'],1024,.3)):
            with tempfile.TemporaryDirectory(prefix='wh2-k5-recovery-cap-') as d:
                paths = [Path(d)/'raw',Path(d)/'error']
                with patch.object(C,'RAW_CAP',cap),patch.object(C.subprocess,'Popen',side_effect=lambda args,**kw:real(helper,**kw)):
                    out,error,code,failure = C.capture('/never-a-codec','0'*64,time.monotonic()+deadline,paths)
                self.assertIsNotNone(failure); self.assertIsNotNone(code); self.assertEqual(error,b'')
                self.assertEqual(out,b'x'*32 if cap==32 else b'prefix\n'); self.assertEqual(paths[0].read_bytes(),out)
                if cap==1024: self.assertEqual(code,-9)

    def test_cross_backend_changed_control_equations_rejected(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k5-recovery-parity-') as d:
            receipt = Path(d)/'receipt'
            receipt.write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL,executables={m:m for m in C.MODES})))
            def capture(executable,claim,deadline,spools):
                records = copy.deepcopy(self.records)
                if executable=='scalar':
                    r = records[1]
                    for a in (2,5): r['arms'][a] = arm_result(r['B'],r['ids'],a,True)
                data = raw(records,executable,claim)
                C.A.publish(spools[0],data); C.A.publish(spools[1],b'')
                return data,b'',0,None
            output = Path(d)/'result'
            with patch.object(C,'OUTPUT',output),patch.object(C,'current'),patch.object(C,'capture',side_effect=capture):
                C.run(receipt)
            analysis = json.loads((output/'analysis.json').read_bytes())
            self.assertEqual(analysis['outcome'],'INVALID')
            self.assertIn('cross-backend complete byte/status parity',analysis['failure'])
            self.assertFalse((output/'asan.raw.jsonl').exists())
            self.assertTrue((output/'COMPLETE.json').exists())


if __name__ == '__main__':
    unittest.main()
