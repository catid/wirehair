"""Synthetic instrument tests; no scientific codec workload or new trace selection."""
import copy
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import time
import unittest
from unittest.mock import patch

import Wh2K8PublicRecoveryR0 as C


def arm_result(e,arm):
    width,tail,ids=e['B'],e['tail'],e['ids']
    rows=tuple(C.coefficient(i) for i in ids)
    first=C.first_success(rows)
    size=7*width+tail
    if arm in (0,3):
        profile=struct.pack('<4sHHQQII',b'WHV2',1,32,0x4b295bbb47f4f9c9,size,width,0)
    elif arm in (1,4):
        profile=struct.pack('<4sHHQQII',b'WHV2',1,32,0x7a9276b85c730ae0,size,width,0)
    else: profile=bytes(32)
    return dict(profile=profile.hex(),rows=bytes(v for row in rows for v in row).hex(),
                packets=[C.packet_hash(width,tail,row,tail if i==7 else width) for i,row in zip(ids,rows)],
                feed=[1]*(first-1)+[0] if first else [1]*len(ids),first=first,
                recoveries=2 if first else 0,counts=[9,9*len(ids),1,first or len(ids),2 if first else 0,10],checked=True)


def synthetic(neutral=True):
    records=[]
    for e in C.roster(neutral):
        arms=[arm_result(e,a) for a in range(3)]
        records.append(dict(type='record',arms=arms+copy.deepcopy(arms),
                            **{k:e[k] for k in ('group','index','B','tail','ids')}))
    return records


def raw(records,mode='native',claim='0'*64,neutral=True):
    header=dict(type='header',protocol=C.PROTOCOL,claim=claim,backend=mode,
                scope='neutral' if neutral else 'retained',retained_raw_sha256=C.D.RAW_SHA,
                features=[0]*4 if mode=='scalar' else [1]*4,
                sources=[C.A.sha(C.message(b,b)) for b in C.WIDTHS])
    return b''.join(C.A.canonical(r) for r in [header]+records+
                    [dict(type='footer',records=48 if neutral else 6260,checked=True)])


class RecoveryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records=synthetic()
        cls.raw=raw(cls.records)

    def reject(self,change,pattern='.'):
        records=copy.deepcopy(self.records);change(records)
        with self.assertRaisesRegex(ValueError,pattern): C.verify(raw(records),'0'*64,'native',True)

    def test_neutral_three_backends_and_tail_shapes(self):
        for mode in C.MODES:
            summary,records=C.verify(raw(self.records,mode),'0'*64,mode,True)
            self.assertEqual(summary,dict(outcome='PASS',neutral_cases=48))
            self.assertEqual(records,self.records)
        self.assertEqual(len(C.roster(True)[:24]),24)
        self.assertEqual([r['index'] for r in C.roster(True)[24:]],list(range(12))+list(range(6132,6144)))
        self.assertEqual({(e['B'],e['tail']) for e in C.roster(True)},
                         {(b,t) for b in C.WIDTHS for t in (b,1)})

    def test_exact_origins_and_horizons(self):
        expected=C.roster()
        self.assertEqual(len(expected),6260)
        self.assertEqual(sum(len(e['ids']) for e in expected),74946)
        for e,origin in zip(expected[6216:],C.retained()[0]['inputs']['origins']):
            self.assertEqual((e['B'],e['tail'],e['ids'],e['origin']),
                             (origin['b'],origin['tail'],origin['ids'],origin))

    def test_every_sealed_coefficient_and_lambda(self):
        for row in C.retained()[1]['rows']:
            self.assertEqual(C.coefficient(row['id']),tuple(row['row']))
        for i in range(8): self.assertEqual(C.coefficient(i),tuple(int(j==i) for j in range(8)))
        self.assertEqual(C.coefficient(8),(96,19,186,153,85,252,7,255))

    def test_tail_packet_hash_and_padding(self):
        for b in C.WIDTHS:
            source=C.message(b,1)
            self.assertEqual(C.packet_hash(b,1,C.coefficient(7),1),C.A.sha(source[7*b:]))
            padded=source+bytes(b-1)
            row=C.coefficient(8)
            expected=bytearray(b)
            for j in range(b):
                for k in range(8): expected[j]^=C.D.R.multiply(row[k],padded[k*b+j])
            self.assertEqual(C.packet_hash(b,1,row,b),C.A.sha(expected))

    def test_shape_chronology_and_tail_rejection(self):
        for change in (lambda r:r.pop(),lambda r:r[0]['arms'].pop(),
                       lambda r:r[0].__setitem__('index',1),lambda r:r[4].__setitem__('tail',2),
                       lambda r:r[-1]['ids'].append(0)):
            self.reject(change)

    def test_all_arm_descriptors_packets_and_ledgers(self):
        for a in range(6):
            self.reject(lambda r:r[0]['arms'][a]['packets'].__setitem__(0,'00'*32),'packet hashes')
            self.reject(lambda r:r[0]['arms'][a]['counts'].__setitem__(0,8),'attempted API ledger')
            self.reject(lambda r:r[0]['arms'][a].__setitem__('rows','ff'+r[0]['arms'][a]['rows'][2:]))
        self.reject(lambda r:r[0]['arms'][0].__setitem__('profile',r[0]['arms'][1]['profile']),'ordinary certified')
        self.reject(lambda r:r[4]['arms'][1].__setitem__('profile',r[0]['arms'][1]['profile']),'WHV2 K8 identity')

    def test_every_arm_rank_status_and_guards(self):
        for a in range(6):
            self.reject(lambda r:r[0]['arms'][a].__setitem__('first',9),'first-success rank')
            self.reject(lambda r:r[0]['arms'][a]['feed'].__setitem__(0,0),'prefix feed statuses')
            self.reject(lambda r:r[0]['arms'][a].__setitem__('recoveries',1),'two byte-exact recoveries')
            self.reject(lambda r:r[0]['arms'][a].__setitem__('checked',False),'guard validation')
        for value in (True,1,7,13): self.reject(lambda r:r[0]['arms'][0].__setitem__('first',value))

    def test_hex_header_truncation_and_scope(self):
        for value in ('00','GG'*32,True,None):
            self.reject(lambda r:r[0]['arms'][0]['packets'].__setitem__(0,value),'hex bytes')
        for value in (self.raw[:-1],self.raw[:100],self.raw.replace(C.D.RAW_SHA.encode(),b'0'*64,1)):
            with self.assertRaises(ValueError): C.verify(value,'0'*64,'native',True)
        with self.assertRaises(ValueError): C.verify(self.raw,'1'*64,'native',True)
        with self.assertRaises(ValueError): C.verify(self.raw,'0'*64,'scalar',True)
        with self.assertRaises(ValueError): C.verify(self.raw,'0'*64,'native',False)

    def test_spent_namespace_does_not_launch(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k8-public-recovery-spent-') as d:
            path=Path(d)/'receipt';path.write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL)))
            with patch.object(C,'OUTPUT',Path(d)),patch.object(C,'current'),patch.object(C,'capture') as capture:
                with self.assertRaises(FileExistsError): C.run(path)
                capture.assert_not_called()

    def test_receipt_cannot_drop_pins_or_rebind_backend(self):
        files,manifests,executables={},{},{}
        def add(path,content):
            pin=dict(path=path,bytes=len(content),sha256=C.A.sha(content))
            files[path]=pin;return pin
        source=add('/synthetic/source.cpp',b'source')
        for mode in C.MODES:
            executable='/synthetic/'+mode+'/recovery_worker';executables[mode]=executable
            artifact=add(executable,mode.encode())
            fixture=add('/synthetic/'+mode+'/fixtures.jsonl',b'fixture')
            path='/synthetic/'+mode+'/manifest.json'
            manifests[path]=C.A.canonical(dict(protocol=C.PROTOCOL,mode=mode,
                                              inputs=[source],artifacts=[artifact,fixture]))
            add(path,manifests[path])
        receipt=dict(protocol=C.PROTOCOL,head='head',executables=executables,
                     environment=dict({k:None for k in C.U.ENV_KEYS},**C.SANITIZERS),
                     pins=sorted(files.values(),key=lambda p:p['path']))
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
                changed=copy.deepcopy(receipt);mutate(changed)
                with self.assertRaises(ValueError):C.current(changed)

    def test_complete_synthetic_nonwinner_stays_fail(self):
        records=synthetic(False)
        summary,observed=C.verify(raw(records,neutral=False),'0'*64,'native')
        self.assertEqual(summary['outcome'],'FAIL')
        self.assertTrue(summary['execution_valid']);self.assertTrue(summary['retained_target_qualified'])
        self.assertFalse(summary['retained_paired_oh0_superior'])
        self.assertEqual(observed,records)
        self.assertEqual(summary['groups'][0]['failure_counts_oh0_to_4'][1],[12,0,0,0,0])

    def test_backend_mismatch_seals_invalid_and_preserves_prefix(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k8-public-backend-parity-') as d:
            bundle=Path(d)/'bundle';receipt=Path(d)/'receipt.json'
            receipt.write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL,
                               executables={m:m for m in C.MODES})))
            calls=[]
            def capture(executable,claim,deadline,spools):
                calls.append(executable)
                output=(executable+'\n').encode()
                C.A.publish(spools[0],output);C.A.publish(spools[1],b'')
                return output,b'',0,None
            def verify(output,claim,mode):
                return dict(outcome='PASS',execution_valid=True),[dict(control_first=8 if mode=='native' else 9)]
            with patch.object(C,'OUTPUT',bundle),patch.object(C,'current'), \
                    patch.object(C,'capture',side_effect=capture),patch.object(C,'verify',side_effect=verify), \
                    patch('builtins.print'):
                C.run(receipt)
                self.assertEqual(calls,['native','scalar'])
                analysis=C.replay()
                self.assertEqual(analysis['outcome'],'INVALID');self.assertFalse(analysis['execution_valid'])
                self.assertIn('cross-backend',analysis['failure'])
                self.assertEqual((bundle/'native.raw.jsonl').read_bytes(),b'native\n')
                self.assertEqual((bundle/'scalar.raw.jsonl').read_bytes(),b'scalar\n')
                self.assertFalse((bundle/'asan.raw.jsonl').exists())
                os.chmod(bundle/'native.raw.jsonl',0o600)
                (bundle/'native.raw.jsonl').write_bytes(b'tampered\n')
                with self.assertRaisesRegex(ValueError,'sealed member identity'):C.replay()

    def test_partial_spool_preserves_existing_file(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k8-public-recovery-spool-') as d:
            paths=[Path(d)/'raw',Path(d)/'error'];paths[1].write_bytes(b'preserved')
            out,err,code,failure=C.capture('/never-run','0'*64,time.monotonic()+1,paths)
            self.assertEqual((out,err,code),(b'',b'',None));self.assertIsNotNone(failure)
            self.assertEqual(paths[0].stat().st_mode&0o777,0o400)
            self.assertEqual(paths[1].read_bytes(),b'preserved')

    def test_output_cap_timeout_and_cleanup(self):
        real=subprocess.Popen
        for code,cap,deadline in [('print("x"*100,flush=True)',32,2),
                                 ('import time;print("prefix",flush=True);time.sleep(2)',1024,.3)]:
            with tempfile.TemporaryDirectory(prefix='wh2-k8-public-recovery-cap-') as d:
                paths=[Path(d)/'raw',Path(d)/'error']
                with patch.object(C,'RAW_CAP',cap),patch.object(C.subprocess,'Popen',
                        side_effect=lambda args,**kw:real(['/usr/bin/python3','-c',code],**kw)):
                    out,err,status,failure=C.capture('/never-a-codec','0'*64,time.monotonic()+deadline,paths)
                self.assertIsNotNone(failure);self.assertIsNotNone(status);self.assertEqual(err,b'')
                self.assertEqual(out,b'x'*32 if cap==32 else b'prefix\n')
                self.assertEqual(paths[0].read_bytes(),out)
                if cap==1024:self.assertEqual(status,-9)


    def test_candidate_descriptor_identity_seed_and_reserved_bytes(self):
        for offset,value in ((0,ord('X')),(8,0),(28,1),(29,1)):
            def mutate(records):
                p=bytearray.fromhex(records[0]['arms'][1]['profile']);p[offset]=value
                records[0]['arms'][1]['profile']=p.hex()
            self.reject(mutate,'installed WHV2 K8 identity')
        def prototype(records):
            records[0]['arms'][1]['profile']=struct.pack('<4sHHQQII',b'WHK8',1,32,
                0x5748324b38544d31,16,2,0).hex()
        self.reject(prototype,'installed WHV2 K8 identity')

    def test_consistently_forged_repeated_control_equation_is_rejected(self):
        def mutate(records):
            record=records[1] # Full B2 low-repair IDs repeat IDs8..11 of record0.
            for arm in (0,3):
                r=record['arms'][arm]
                encoded=bytearray.fromhex(r['rows']);encoded[0]^=1
                rows=tuple(tuple(encoded[i:i+8]) for i in range(0,len(encoded),8))
                r['rows']=encoded.hex()
                r['packets']=[C.packet_hash(record['B'],record['tail'],row,
                    record['tail'] if pid==7 else record['B']) for pid,row in zip(record['ids'],rows)]
                first=C.first_success(rows)
                r.update(first=first,feed=[1]*(first-1)+[0] if first else [1]*len(rows),
                         recoveries=2 if first else 0,
                         counts=[9,9*len(rows),1,first or len(rows),2 if first else 0,10])
        self.reject(mutate,'repeated descriptor/ID equation and payload')

    def test_prototype_parity_never_rewrites_reported_descriptors(self):
        # Small mocked archive; the real full roster count remains enforced.
        records=copy.deepcopy(self.records)
        records=(records*131)[:6260]
        old=copy.deepcopy(records)
        for r in old:
            for a in (1,4):
                r['arms'][a]['profile']=struct.pack('<4sHHQQII',b'WHK8',1,32,
                    0x5748324b38544d31,7*r['B']+r['tail'],r['B'],0).hex()
        old_raw=b''.join(C.A.canonical(r) for r in [{}]+old+[{}])
        raw_path=C.PROTOTYPE/'native.raw.jsonl'
        pin=dict(path=str(raw_path),bytes=len(old_raw),sha256=C.A.sha(old_raw))
        complete=C.A.canonical(dict(protocol='wirehair.wh2.k8-serialized-recovery-r0',
                                   outcome='PASS',files=[pin]))
        before=copy.deepcopy(records)
        with patch.object(C,'PROTOTYPE_COMPLETE_SHA',C.A.sha(complete)), \
                patch.object(C,'PROTOTYPE_RAW_SHA',C.A.sha(old_raw)), \
                patch.object(C.A,'read_regular',side_effect=lambda p,cap:old_raw if p==raw_path else complete), \
                patch.object(C.U,'pin',return_value=pin):
            C.retained_parity(records)
            self.assertEqual(records,before)
            for arm in range(6):
                changed=copy.deepcopy(records);changed[0]['arms'][arm]['checked']=False
                with self.assertRaisesRegex(ValueError,'parity'):C.retained_parity(changed)


class DecisionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.expected=C.roster()
        cls.records=[dict(arms=[dict(first=8) for _ in range(6)]) for _ in cls.expected]
        # A purely synthetic control deficiency, not a native codec observation.
        for a in (0,2,3,5):cls.records[0]['arms'][a]['first']=9

    def result(self,change=None):
        rows=copy.deepcopy(self.records)
        if change:change(rows)
        return C.summarize(rows,self.expected)

    def test_three_separate_predicates_and_claim_limits(self):
        s=self.result();self.assertEqual(s['outcome'],'PASS')
        for key in ('execution_valid','retained_target_qualified','retained_paired_oh0_superior'):self.assertTrue(s[key])
        for key in ('fresh_sample','speed_claimed','all_K_claimed','production_promotion_claimed'):self.assertFalse(s[key])
        self.assertEqual([g['cases'] for g in s['groups']],[6144,72,44])

    def test_valid_nonwinning_is_fail_not_invalid_or_pass(self):
        def tie(rows):
            for a in (0,2,3,5):rows[0]['arms'][a]['first']=8
        s=self.result(tie)
        self.assertEqual(s['outcome'],'FAIL');self.assertTrue(s['execution_valid'])
        self.assertTrue(s['retained_target_qualified']);self.assertFalse(s['retained_paired_oh0_superior'])

    def test_each_pair_is_required(self):
        for a in (0,2,3,5):
            s=self.result(lambda rows:rows[0]['arms'][a].__setitem__('first',8))
            self.assertEqual(s['outcome'],'FAIL')

    def test_cell_integer_threshold_and_hard_history_requirements(self):
        indices=[i for i,e in enumerate(self.expected) if e['group']==0 and e['B']==2 and
                 e['schedule']==self.expected[0]['schedule']][:6]
        for count,qualified in ((5,True),(6,False)):
            def change(rows):
                for i in indices[:count]:rows[i]['arms'][1]['first']=9
            self.assertEqual(self.result(change)['retained_target_qualified'],qualified)
        for index in (6144,6216):
            for a in (1,4):
                s=self.result(lambda rows:rows[index]['arms'][a].__setitem__('first',0))
                self.assertFalse(s['retained_target_qualified']);self.assertEqual(s['outcome'],'FAIL')

    def test_introductions_and_history_horizons_are_preserved(self):
        def change(rows):
            rows[1]['arms'][1]['first']=9
            rows[6216]['arms'][1]['first']=0
        s=self.result(change)
        key=[self.expected[1][k] for k in ('index','B','tail')]
        self.assertIn(key,s['groups'][0]['paired'][0]['introduced'][0])
        for group in s['groups'][:2]+s['cells']:
            counts=group['failure_counts_oh0_to_4']
            for p,(a,b) in zip(group['paired'],C.PAIRED):
                for oh in range(5):self.assertEqual(len(p['fixed'][oh])-len(p['introduced'][oh]),counts[b][oh]-counts[a][oh])
        history=s['groups'][2]['paired'][0]
        length=len(self.expected[6216]['ids'])
        for oh in range(5):self.assertEqual(bool(history['introduced'][oh]),8+oh<=length)


if __name__=='__main__':unittest.main()
