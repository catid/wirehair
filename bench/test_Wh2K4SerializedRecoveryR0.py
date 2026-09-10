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

import Wh2K4SerializedRecoveryR0 as C


def arm_result(e,arm):
    width,tail,ids=e['B'],e['tail'],e['ids']
    rows=tuple(C.coefficient(i) for i in ids)
    first=C.first_success(rows)
    size=3*width+tail
    if arm in (0,3):
        profile=struct.pack('<4sHHQQII',b'WHV2',1,32,0x4b295bbb47f4f9c9,size,width,0)
    elif arm in (1,4):
        profile=struct.pack('<4sHHQQII',b'WHK4',1,32,0x5748324b34544d31,size,width,0)
    else: profile=bytes(32)
    return dict(profile=profile.hex(),rows=bytes(v for row in rows for v in row).hex(),
                packets=[C.packet_hash(width,tail,row,tail if i==3 else width) for i,row in zip(ids,rows)],
                feed=[1]*(first-1)+[0] if first else [1]*len(ids),first=first,
                recoveries=2 if first else 0,counts=[5,5*len(ids),1,first or len(ids),2 if first else 0,6],checked=True)


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
                    [dict(type='footer',records=len(records) if neutral else 6254,checked=True)])


def replace_rows(record, arm, rows):
    """Coherent forged control output: recompute hashes, endpoints and ledgers."""
    e = record; rows = tuple(tuple(r) for r in rows)
    r = e['arms'][arm]; first = C.first_success(rows)
    r.update(rows=bytes(v for row in rows for v in row).hex(),
             packets=[C.packet_hash(e['B'],e['tail'],row,e['tail'] if i==3 else e['B'])
                      for i,row in zip(e['ids'],rows)],
             first=first,feed=[1]*(first-1)+[0] if first else [1]*len(rows),
             recoveries=2 if first else 0,counts=[5,5*len(rows),1,first or len(rows),2 if first else 0,6])
    e['arms'][arm+3] = copy.deepcopy(r)


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
        self.assertEqual(len(expected),6254)
        self.assertEqual(sum(len(e['ids']) for e in expected),49891)
        for e,origin in zip(expected[6216:],C.retained()[0]['inputs']['origins']):
            self.assertEqual((e['B'],e['tail'],e['ids'],e['origin']),
                             (origin['b'],origin['tail'],origin['ids'],origin))

    def test_every_sealed_coefficient_and_lambda(self):
        for row in C.retained()[1]['rows']:
            self.assertEqual(C.coefficient(row['id']),tuple(row['row']))
        for i in range(4): self.assertEqual(C.coefficient(i),tuple(int(j==i) for j in range(4)))
        self.assertEqual(C.coefficient(4),(64,120,54,15))

    def test_tail_packet_hash_and_padding(self):
        for b in C.WIDTHS:
            source=C.message(b,1)
            self.assertEqual(C.packet_hash(b,1,C.coefficient(3),1),C.A.sha(source[3*b:]))
            padded=source+bytes(b-1)
            row=C.coefficient(4)
            expected=bytearray(b)
            for j in range(b):
                for k in range(4): expected[j]^=C.D.R.multiply(row[k],padded[k*b+j])
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
        self.reject(lambda r:r[4]['arms'][1].__setitem__('profile',r[0]['arms'][1]['profile']),'WHK4 identity')

    def test_every_arm_rank_status_and_guards(self):
        for a in range(6):
            self.reject(lambda r:r[0]['arms'][a].__setitem__('first',5),'first-success rank')
            self.reject(lambda r:r[0]['arms'][a]['feed'].__setitem__(0,0),'prefix feed statuses')
            self.reject(lambda r:r[0]['arms'][a].__setitem__('recoveries',1),'two byte-exact recoveries')
            self.reject(lambda r:r[0]['arms'][a].__setitem__('checked',False),'guard validation')
        for value in (True,1,3,9): self.reject(lambda r:r[0]['arms'][0].__setitem__('first',value))

    def test_coherent_repeated_control_rows_are_rejected(self):
        for arm in (0,2):
            for record_index,row_index in ((3,1),(1,0)):
                def mutate(records):
                    record = records[record_index]
                    rows = [list(C.coefficient(i)) for i in record['ids']]
                    rows[row_index][1] ^= 1
                    replace_rows(record,arm,rows)
                self.reject(mutate,'repeated packet row identity')

    def test_control_endpoints_may_be_before_or_after_candidate(self):
        deficient = next(e['ids'] for e in C.roster()[:6144]
                         if len(set(e['ids']))==8 and C.candidate_first(tuple(e['ids']))==5)
        for ids,expected in ((list(range(8)),[5,4]),(deficient,[4,5])):
            e = dict(group=3,index=0,B=2,tail=2,ids=ids)
            with patch.object(C,'roster',return_value=[e]):
                records = synthetic()
                # Independent control equations are synthetic, never codec evidence.
                mapped = [0,1,2,2,3,4,5,6] if expected[0]==5 else list(range(8))
                controls = [C.coefficient(i) for i in mapped]
                replace_rows(records[0],0,controls)
                _, observed = C.verify(raw(records),'0'*64,'native',True)
                self.assertEqual([observed[0]['arms'][a]['first'] for a in (0,1)], expected)

    def test_hex_header_truncation_and_scope(self):
        for value in ('00','GG'*32,True,None):
            self.reject(lambda r:r[0]['arms'][0]['packets'].__setitem__(0,value),'hex bytes')
        for value in (self.raw[:-1],self.raw[:100],self.raw.replace(C.D.RAW_SHA.encode(),b'0'*64,1)):
            with self.assertRaises(ValueError): C.verify(value,'0'*64,'native',True)
        with self.assertRaises(ValueError): C.verify(self.raw,'1'*64,'native',True)
        with self.assertRaises(ValueError): C.verify(self.raw,'0'*64,'scalar',True)
        with self.assertRaises(ValueError): C.verify(self.raw,'0'*64,'native',False)

    def test_spent_namespace_does_not_launch(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k4-recovery-spent-') as d:
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
            # Coherently omit the consumed fixture from BOTH metadata layers.
            fixture='/synthetic/native/fixtures.jsonl';files.pop(fixture)
            path='/synthetic/native/manifest.json';manifest=C.A.decode(manifests[path])
            manifest['artifacts']=[p for p in manifest['artifacts'] if p['path']!=fixture]
            manifests[path]=C.A.canonical(manifest);add(path,manifests[path])
            changed=copy.deepcopy(receipt);changed['pins']=sorted(files.values(),key=lambda p:p['path'])
            with self.assertRaisesRegex(ValueError,'pinned neutral fixture artifact'):C.current(changed)

    def test_complete_synthetic_nonwinner_stays_fail(self):
        records=synthetic(False)
        summary,observed=C.verify(raw(records,neutral=False),'0'*64,'native')
        self.assertEqual(summary['outcome'],'FAIL')
        self.assertTrue(summary['execution_valid']);self.assertTrue(summary['retained_target_qualified'])
        self.assertFalse(summary['retained_paired_oh0_superior'])
        self.assertEqual(observed,records)
        self.assertEqual(summary['groups'][0]['failure_counts_oh0_to_4'][1],[9,0,0,0,0])

    def test_backend_mismatch_seals_invalid_and_preserves_prefix(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k4-backend-parity-') as d:
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
                return dict(outcome='PASS',execution_valid=True),[dict(control_first=4 if mode=='native' else 5)]
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
        with tempfile.TemporaryDirectory(prefix='wh2-k4-recovery-spool-') as d:
            paths=[Path(d)/'raw',Path(d)/'error'];paths[1].write_bytes(b'preserved')
            out,err,code,failure=C.capture('/never-run','0'*64,time.monotonic()+1,paths)
            self.assertEqual((out,err,code),(b'',b'',None));self.assertIsNotNone(failure)
            self.assertEqual(paths[0].stat().st_mode&0o777,0o400)
            self.assertEqual(paths[1].read_bytes(),b'preserved')

    def test_worker_environment_and_cwd_are_explicit(self):
        real = subprocess.Popen
        hostile = dict(LD_AUDIT='/bad/audit',LD_DEBUG='all',LD_DEBUG_OUTPUT='/bad/log',
                       CPATH='/bad/includes',PYTHONPATH='/bad/python',HOME='/bad/home',
                       **C.SANITIZERS)
        calls = []
        def spawn(argv,**kw):
            calls.append(kw)
            return real(['/usr/bin/python3','-c','print("neutral")'],**kw)
        with tempfile.TemporaryDirectory(prefix='wh2-k4-env-') as d, \
                patch.dict(os.environ,hostile,clear=True),patch.object(C.subprocess,'Popen',side_effect=spawn):
            out,err,code,failure=C.capture('/not-a-codec','0'*64,time.monotonic()+5,[Path(d)/'raw',Path(d)/'err'])
        self.assertEqual((out,err,code,failure),(b'neutral\n',b'',0,None))
        self.assertEqual(calls[0]['env'],dict(PATH='/usr/bin:/bin',LANG='C',LC_ALL='C',TZ='UTC',**C.SANITIZERS))
        self.assertEqual(calls[0]['cwd'],C.ROOT)

    def test_expired_capture_and_preflight_never_claim_or_launch(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k4-expired-') as d,patch.object(C.subprocess,'Popen') as popen:
            paths=[Path(d)/'raw',Path(d)/'err']
            out,err,status,failure=C.capture('/not-a-codec','0'*64,time.monotonic()-1,paths)
            self.assertEqual((out,err,status),(b'',b'',None));self.assertIn('deadline',failure)
            self.assertTrue(all(not p.exists() for p in paths));popen.assert_not_called()
            receipt=Path(d)/'receipt';receipt.write_bytes(C.A.canonical({}))
            bundle=Path(d)/'bundle'
            with patch.object(C,'OUTPUT',bundle),patch.object(C,'current'), \
                    patch.object(C.time,'monotonic',side_effect=[0,1000]),patch.object(C,'capture') as capture:
                with self.assertRaisesRegex(ValueError,'deadline'):C.run(receipt)
                self.assertFalse(bundle.exists());capture.assert_not_called()

    def test_spool_exact_caps_short_and_failed_writes(self):
        real_spawn,real_write=subprocess.Popen,os.write
        for stream in (0,1):
            for size in (31,32,33):
                code='import os;os.write(%d,b"x"*%d)'%(stream+1,size)
                with tempfile.TemporaryDirectory(prefix='wh2-k4-spool-') as d, \
                        patch.object(C,'RAW_CAP',32),patch.object(C,'ERR_CAP',32), \
                        patch.object(C.subprocess,'Popen',side_effect=lambda argv,**kw:real_spawn(['/usr/bin/python3','-c',code],**kw)), \
                        patch.object(C.os,'write',side_effect=lambda fd,data:real_write(fd,data[:3])):
                    paths=[Path(d)/'raw',Path(d)/'err']
                    out,err,status,failure=C.capture('/not-a-codec','0'*64,time.monotonic()+5,paths)
                    self.assertEqual((out,err)[stream],b'x'*min(size,32))
                    self.assertEqual(paths[stream].read_bytes(),(out,err)[stream])
                    self.assertEqual(failure is None,size<=32)
        for action in ('write','fsync','close'):
            real_action=getattr(os,action);started=[False]
            def spawn(argv,**kw):
                child=real_spawn(['/usr/bin/python3','-c','print("x")'],**kw)
                started[0]=True
                return child
            def fail(*args):
                if not started[0]:return real_action(*args)
                if action=='write':return 0
                if action=='close':real_action(*args)
                raise OSError('injected '+action)
            with tempfile.TemporaryDirectory(prefix='wh2-k4-spool-fault-') as d:
                with patch.object(C.subprocess,'Popen',side_effect=spawn), \
                        patch.object(C.os,action,side_effect=fail):
                    _,_,_,failure=C.capture('/not-a-codec','0'*64,time.monotonic()+5,[Path(d)/'raw',Path(d)/'err'])
                self.assertIsNotNone(failure)

    def test_output_cap_timeout_and_cleanup(self):
        real=subprocess.Popen
        for code,cap,deadline in [('print("x"*100,flush=True)',32,2),
                                 ('import time;print("prefix",flush=True);time.sleep(2)',1024,.3)]:
            with tempfile.TemporaryDirectory(prefix='wh2-k4-recovery-cap-') as d:
                paths=[Path(d)/'raw',Path(d)/'error']
                with patch.object(C,'RAW_CAP',cap),patch.object(C.subprocess,'Popen',
                        side_effect=lambda args,**kw:real(['/usr/bin/python3','-c',code],**kw)):
                    out,err,status,failure=C.capture('/never-a-codec','0'*64,time.monotonic()+deadline,paths)
                self.assertIsNotNone(failure);self.assertIsNotNone(status);self.assertEqual(err,b'')
                self.assertEqual(out,b'x'*32 if cap==32 else b'prefix\n')
                self.assertEqual(paths[0].read_bytes(),out)
                if cap==1024:self.assertEqual(status,-9)


class DecisionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.expected=C.roster()
        cls.records=[dict(arms=[dict(first=4) for _ in range(6)]) for _ in cls.expected]
        # A purely synthetic control deficiency, not a native codec observation.
        for a in (0,2,3,5):cls.records[0]['arms'][a]['first']=5

    def result(self,change=None):
        rows=copy.deepcopy(self.records)
        if change:change(rows)
        return C.summarize(rows,self.expected)

    def test_three_separate_predicates_and_claim_limits(self):
        s=self.result();self.assertEqual(s['outcome'],'PASS')
        for key in ('execution_valid','retained_target_qualified','retained_paired_oh0_superior'):self.assertTrue(s[key])
        for key in ('fresh_sample','speed_claimed','all_K_claimed','production_promotion_claimed'):self.assertFalse(s[key])
        self.assertEqual([g['cases'] for g in s['groups']],[6144,72,38])

    def test_valid_nonwinning_is_fail_not_invalid_or_pass(self):
        def tie(rows):
            for a in (0,2,3,5):rows[0]['arms'][a]['first']=4
        s=self.result(tie)
        self.assertEqual(s['outcome'],'FAIL');self.assertTrue(s['execution_valid'])
        self.assertTrue(s['retained_target_qualified']);self.assertFalse(s['retained_paired_oh0_superior'])

    def test_each_pair_is_required(self):
        for a in (0,2,3,5):
            s=self.result(lambda rows:rows[0]['arms'][a].__setitem__('first',4))
            self.assertEqual(s['outcome'],'FAIL')

    def test_cell_integer_threshold_and_hard_history_requirements(self):
        indices=[i for i,e in enumerate(self.expected) if e['group']==0 and e['B']==2 and
                 e['schedule']==self.expected[0]['schedule']][:6]
        for a in (1,4):
            for count,qualified in ((5,True),(6,False)):
                def change(rows):
                    for i in indices[:count]:rows[i]['arms'][a]['first']=5
                self.assertEqual(self.result(change)['retained_target_qualified'],qualified)
        for index in (6144,6216):
            for a in (1,4):
                s=self.result(lambda rows:rows[index]['arms'][a].__setitem__('first',0))
                self.assertFalse(s['retained_target_qualified']);self.assertEqual(s['outcome'],'FAIL')

    def test_introductions_and_history_horizons_are_preserved(self):
        def change(rows):
            rows[1]['arms'][1]['first']=5
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
        for oh in range(5):self.assertEqual(bool(history['introduced'][oh]),4+oh<=length)

    def test_historical_counts_never_invent_horizons(self):
        def change(rows):
            for r,e in zip(rows,self.expected):
                if e['group']==2:
                    for arm in r['arms']:arm['first']=0
        summary=self.result(change);history=summary['groups'][2]
        self.assertEqual(history['eligible_cases_oh0_to_4'],[38,9,2,0,0])
        self.assertEqual(history['failure_counts_oh0_to_4'],[[38,9,2,0,0]]*6)
        self.assertEqual(history['unresolved'],[38]*6)

    def test_two_fixes_one_introduction_is_a_valid_net_win(self):
        def change(rows):
            for a in (0,2,3,5):rows[1]['arms'][a]['first']=5
            for a in (1,4):rows[2]['arms'][a]['first']=5
        summary=self.result(change)
        self.assertEqual(summary['outcome'],'PASS')
        for p in summary['groups'][0]['paired']:
            self.assertEqual((len(p['fixed'][0]),len(p['introduced'][0])),(2,1))


class InvalidReplayTest(unittest.TestCase):
    def bundle(self,d,analysis=None,alter_bytes=None,alter_members=None):
        root=Path(d)
        stored=dict(outcome='INVALID',execution_valid=False,failure='injected',
                    protocol=C.PROTOCOL,elapsed_seconds=1.0)
        if analysis:analysis(stored)
        encoded=C.A.canonical(stored)
        if alter_bytes:encoded=alter_bytes(encoded)
        (root/'CLAIM.json').write_bytes(C.A.canonical({}))
        (root/'analysis.json').write_bytes(encoded)
        (root/'native.raw.jsonl').write_bytes(b'prefix\n')
        (root/'native.stderr.txt').write_bytes(b'failed\n')
        members=[C.U.pin(p) for p in sorted(root.iterdir())]
        if alter_members:alter_members(members)
        (root/'COMPLETE.json').write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL,outcome='INVALID',files=members)))
        return root

    def test_strict_invalid_schema_and_finite_timing(self):
        changes=[(None,None,None),
                 (lambda a:a.__setitem__('execution_valid',True),None,'INVALID execution'),
                 (lambda a:a.__setitem__('failure',False),None,'INVALID failure'),
                 (lambda a:a.__setitem__('retained_target_qualified',True),None,'INVALID analysis schema')]
        changes += [(lambda a,v=v:a.__setitem__('elapsed_seconds',v),None,'elapsed time') for v in (-1,True,'1')]
        changes.append((None,lambda b:b.replace(b'"elapsed_seconds":1.0',b'"elapsed_seconds":1e999'),'elapsed time'))
        for change,raw_change,error in changes:
            with tempfile.TemporaryDirectory(prefix='wh2-k4-invalid-replay-') as d:
                root=self.bundle(d,change,raw_change)
                with patch.object(C,'OUTPUT',root),patch.object(C,'current'),patch.object(C,'verify') as verify:
                    if error:
                        with self.assertRaisesRegex(ValueError,error):C.replay()
                    else:self.assertEqual(C.replay()['outcome'],'INVALID')
                    verify.assert_not_called()

    def test_invalid_replay_member_and_bundle_limits(self):
        for name,cap in (('native.raw.jsonl',C.RAW_CAP),('native.stderr.txt',C.ERR_CAP)):
            def change(members):
                next(p for p in members if Path(p['path']).name==name)['bytes']=cap+1
            with tempfile.TemporaryDirectory(prefix='wh2-k4-invalid-cap-') as d:
                root=self.bundle(d,alter_members=change)
                with patch.object(C,'OUTPUT',root),patch.object(C,'current'):
                    with self.assertRaisesRegex(ValueError,'sealed member cap'):C.replay()
        with tempfile.TemporaryDirectory(prefix='wh2-k4-invalid-bundle-') as d:
            root=self.bundle(d)
            total=sum(p['bytes'] for p in C.A.decode((root/'COMPLETE.json').read_bytes())['files'])
            with patch.object(C,'OUTPUT',root),patch.object(C,'current'),patch.object(C,'BUNDLE_CAP',65536+total):
                with self.assertRaisesRegex(ValueError,'bundle cap'):C.replay()


if __name__=='__main__':unittest.main()
