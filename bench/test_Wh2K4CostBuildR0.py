"""Synthetic cost build/receipt/capture checks; never runs codec science."""
import copy
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest
from unittest.mock import patch

import Wh2K4SerializedCostR0 as C


class BuildTest(unittest.TestCase):
    def fixture(self):
        return dict(type='header', protocol=C.PROTOCOL, claim='0'*64,
                    batch=128, identity_hex=b'neutral'.hex(), fixtures=[],
                    prelude=dict(clocks=[100,100,102,112,108,120], before=[0]*4, after=[0]*4))

    def fake_receipt(self):
        files, data, manifests, executables = {}, {}, {}, {}
        def add(path, raw):
            path = str(path)
            data[path] = raw
            files[path] = dict(path=path, bytes=len(raw), sha256=C.A.sha(raw))
            return files[path]
        source = add('/synthetic/source.cpp', b'source')
        for mode in C.MODES:
            base = Path('/synthetic')/mode
            artifacts = [add(base/name, C.A.canonical(self.fixture()) if name == 'fixtures.json' else
                             C.A.canonical(dict(producing_source_closure=True)) if name == 'qualified-library.json' else name.encode())
                         for name in ('cost_worker','fixtures.json','contract.json','qualified-library.json')]
            if mode == 'native': artifacts.append(add(base/'target.json', b'target'))
            manifest = dict(protocol=C.PROTOCOL, mode=mode, inputs=[source], artifacts=artifacts)
            manifests[mode] = manifest
            add(base/'manifest.json', C.A.canonical(manifest))
            executables[mode] = str(base/'cost_worker')
        receipt = dict(protocol=C.PROTOCOL, head='head', executable=executables['native'],
                       executables=executables,
                       environment=dict({k: None for k in C.R.ENV_KEYS}, **C.SANITIZERS),
                       pins=sorted(files.values(), key=lambda p:p['path']))
        return receipt, files, data, manifests

    def check_receipt(self, receipt, files, data, producer_prerequisite_mocked=True):
        # These generic schema tests isolate the unrelated producer prerequisite.
        # No real producer is qualified. Direct guard tests below never mock it.
        producer_check=(lambda proof: None) if producer_prerequisite_mocked else C.R.verify_qualified_library
        with patch.dict(os.environ, C.SANITIZERS, clear=True), \
                patch.object(C.R,'verify_qualified_library',side_effect=producer_check), \
                patch.object(C, 'pin', side_effect=lambda p:files[str(p)]), \
                patch.object(C, 'command', return_value=b'head\n'), \
                patch.object(C.A, 'read_regular', side_effect=lambda p,cap:data[str(p)]), \
                patch.object(C, 'prior_header', return_value=self.fixture()), \
                patch.object(C, 'verify_fixtures') as verify:
            C.current(receipt)
            self.assertEqual(verify.call_count, 3)

    def test_receipt_requires_all_backend_closures(self):
        receipt, files, data, _ = self.fake_receipt()
        self.check_receipt(receipt, files, data)
        for mutate in (lambda r:r['pins'].pop(), lambda r:r['pins'].append(r['pins'][0]),
                       lambda r:r['executables'].pop('asan'),
                       lambda r:r.__setitem__('executable',r['executables']['scalar']),
                       lambda r:r['executables'].__setitem__('scalar',r['executables']['native']),
                       lambda r:r.__setitem__('head','other'),
                       lambda r:r.__setitem__('protocol','wrong')):
            changed=copy.deepcopy(receipt); mutate(changed)
            with self.assertRaises(ValueError): self.check_receipt(changed,files,data)

    def test_required_artifact_cannot_be_reclassified_or_omitted(self):
        for mode in C.MODES:
            names = ['cost_worker','fixtures.json','contract.json','qualified-library.json']
            if mode == 'native': names.append('target.json')
            for name in names:
                for keep_as_input in (False, True):
                    receipt, files, data, manifests = self.fake_receipt()
                    path = '/synthetic/'+mode+'/'+name
                    manifest = manifests[mode]
                    artifact = next(p for p in manifest['artifacts'] if p['path']==path)
                    manifest['artifacts'].remove(artifact)
                    if keep_as_input: manifest['inputs'].append(artifact)
                    else:
                        receipt['pins'].remove(artifact)
                        del files[path]
                    manifest_path='/synthetic/'+mode+'/manifest.json'
                    raw=C.A.canonical(manifest); data[manifest_path]=raw
                    files[manifest_path].update(bytes=len(raw),sha256=C.A.sha(raw))
                    with self.subTest(mode=mode,name=name,input=keep_as_input), self.assertRaises(ValueError):
                        self.check_receipt(receipt,files,data)

    def test_artifact_hash_disagreement_rejected(self):
        receipt, files, data, manifests = self.fake_receipt()
        manifest=copy.deepcopy(manifests['native'])
        next(p for p in manifest['artifacts'] if p['path'].endswith('/target.json'))['sha256']='f'*64
        path='/synthetic/native/manifest.json'; raw=C.A.canonical(manifest)
        data[path]=raw; files[path].update(bytes=len(raw),sha256=C.A.sha(raw))
        with self.assertRaisesRegex(ValueError,'artifact identity'):
            self.check_receipt(receipt,files,data)

    def test_fresh_producing_assertion_never_qualifies(self):
        for fresh in (True,False,1,0,None,'qualified'):
            with self.subTest(fresh=fresh), self.assertRaisesRegex(ValueError,'K4 R0'):
                C.R.verify_qualified_library(dict(fresh_neutral=fresh,producing_source_closure=True))
        for closure in (True,False,None,1,'true'):
            with self.subTest(closure=closure), self.assertRaises(ValueError):
                C.R.verify_qualified_library(dict(producing_source_closure=closure))

    def test_receipt_rejects_even_hash_consistent_fresh_proof(self):
        receipt,files,data,manifests=self.fake_receipt()
        proof_path='/synthetic/native/qualified-library.json'
        raw=C.A.canonical(dict(producing_source_closure=True,fresh_neutral=True))
        data[proof_path]=raw; files[proof_path].update(bytes=len(raw),sha256=C.A.sha(raw))
        path='/synthetic/native/manifest.json'
        raw=C.A.canonical(manifests['native']); data[path]=raw
        files[path].update(bytes=len(raw),sha256=C.A.sha(raw))
        with self.assertRaisesRegex(ValueError,'K4 R0'):
            self.check_receipt(receipt,files,data,producer_prerequisite_mocked=False)

    def test_receipt_without_fresh_marker_is_still_blocked(self):
        receipt,files,data,_=self.fake_receipt()
        with self.assertRaisesRegex(ValueError,'K4 R0'):
            self.check_receipt(receipt,files,data,producer_prerequisite_mocked=False)

    def test_historical_recovery_reader_restored_byte_exactly(self):
        original=C.command(['git','show','adb71c2:bench/Wh2K4RecoveryBuildR0.py'])
        self.assertEqual(original,(C.ROOT/'bench/Wh2K4RecoveryBuildR0.py').read_bytes())

    def test_environment_values_are_not_optional(self):
        receipt,files,data,_=self.fake_receipt()
        for key in C.ENV_KEYS:
            changed=copy.deepcopy(receipt); changed['environment'][key]='different'
            with self.assertRaisesRegex(ValueError,'runtime environment'):
                self.check_receipt(changed,files,data)

    def test_expired_capture_never_creates_files_or_child(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k4-cost-expired-') as d:
            paths=[Path(d)/'raw',Path(d)/'error']
            with patch.object(C.subprocess,'Popen') as child:
                raw,err,code,failure=C.capture('/not-run','0'*64,time.monotonic()-1,paths)
                child.assert_not_called()
            self.assertEqual((raw,err,code),(b'',b'',None))
            self.assertIsNotNone(failure)
            self.assertFalse(any(p.exists() for p in paths))

    def test_post_open_deadline_rejects_before_child(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k4-cost-deadline-') as d:
            paths=[Path(d)/'raw',Path(d)/'error']
            with patch.object(C.A,'time_left',side_effect=[1,ValueError('deadline')]), \
                    patch.object(C.subprocess,'Popen') as child:
                result=C.capture('/not-run','0'*64,1,paths)
                child.assert_not_called()
            self.assertIsNotNone(result[3])
            self.assertTrue(all(p.stat().st_mode&0o777==0o400 for p in paths))

    def test_expired_preflight_does_not_spend_namespace(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k4-cost-preflight-') as d:
            base=Path(d); receipt=base/'receipt'
            receipt.write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL)))
            with patch.object(C,'OUTPUT',base/'science'), patch.object(C,'current'), \
                    patch.object(C.time,'monotonic',side_effect=[0,481]), \
                    patch.object(C,'capture') as capture:
                with self.assertRaisesRegex(ValueError,'deadline'): C.run(receipt)
                capture.assert_not_called()
            self.assertFalse((base/'science').exists())

    def invalid_bundle(self, base, change=None, complete_change=None):
        analysis=dict(protocol=C.PROTOCOL,outcome='INVALID',failure='expected',elapsed_seconds=1.0)
        if change: change(analysis)
        # Direct JSON spelling covers finite-parser overflow without emitting NaN.
        import json
        raw=json.dumps(analysis).replace('Infinity','1e999').encode()+b'\n'
        (base/'analysis.json').write_bytes(raw)
        (base/'CLAIM.json').write_bytes(C.A.canonical(dict(protocol=C.PROTOCOL)))
        for path in base.iterdir(): path.chmod(0o400)
        complete=dict(protocol=C.PROTOCOL,outcome='INVALID',files=[C.pin(p) for p in sorted(base.iterdir())])
        if complete_change: complete_change(complete)
        (base/'COMPLETE.json').write_bytes(C.A.canonical(complete)); (base/'COMPLETE.json').chmod(0o400)

    def test_invalid_replay_rejects_claim_and_elapsed_forgery(self):
        mutations=[lambda a:a.update(speed_qualified=True),lambda a:a.update(failure=False)]
        mutations += [lambda a,v=v:a.update(elapsed_seconds=v) for v in (True,-1,float('inf'),'1')]
        for change in mutations:
            with tempfile.TemporaryDirectory(prefix='wh2-k4-cost-invalid-') as d:
                base=Path(d); self.invalid_bundle(base,change=change)
                with patch.object(C,'OUTPUT',base),patch.object(C,'current'),self.assertRaises(ValueError):
                    C.replay()

    def test_invalid_replay_rejects_schema_and_caps_before_member_reads(self):
        changes=[lambda c:c.update(extra=True),
                 lambda c:c['files'][0].update(bytes=1024**2+1)]
        for change in changes:
            with tempfile.TemporaryDirectory(prefix='wh2-k4-cost-seal-') as d:
                base=Path(d);self.invalid_bundle(base,complete_change=change)
                real=C.A.read_regular; read=[]
                def limited(path,cap):
                    read.append(Path(path).name)
                    return real(path,cap)
                with patch.object(C,'OUTPUT',base),patch.object(C,'current'), \
                        patch.object(C.A,'read_regular',side_effect=limited),self.assertRaises(ValueError):
                    C.replay()
                self.assertEqual(read,['COMPLETE.json'])

    def test_output_device_attempt_is_retained_before_rejection(self):
        for sink in ('full','broken-pipe'):
            for response in (subprocess.CompletedProcess([],0,b'',b'wrong'),
                             subprocess.TimeoutExpired(['/not-codec'],60,stderr=b'prefix'),
                             OSError('launch failed')):
                with tempfile.TemporaryDirectory(prefix='wh2-k4-cost-output-') as d:
                    base=Path(d); commands=[]
                    kwargs={'side_effect':response} if isinstance(response,Exception) else {'return_value':response}
                    with patch.object(C.R.subprocess,'run',**kwargs),self.assertRaises(ValueError):
                        C.R.output_device_check('/not-codec',sink,base,0,commands)
                    evidence=C.A.decode((base/'output-device-0.json').read_bytes())
                    self.assertEqual(evidence['sink'],sink)
                    self.assertEqual(evidence['argv'],['/not-codec','--neutral-publication','success'])
                    self.assertEqual(evidence['cwd'],str(C.ROOT))
                    self.assertEqual(len(commands),1)
                    stderr=(base/'output-device-0.stderr').read_bytes()
                    self.assertEqual(stderr,b'prefix' if isinstance(response,subprocess.TimeoutExpired) else
                                     b'' if isinstance(response,OSError) else b'wrong')

    def test_capture_explicit_child_environment_and_cwd(self):
        real=subprocess.Popen
        helper=['/usr/bin/python3','-c',
                'import os; print(os.getcwd()); print(os.environ.get("CPATH", "absent"))']
        observed={}
        def launch(argv,**kwargs):
            observed.update(kwargs)
            return real(helper,**kwargs)
        with tempfile.TemporaryDirectory(prefix='wh2-k4-cost-env-') as d, \
                patch.dict(os.environ,dict(C.SANITIZERS,CPATH='/hostile'),clear=True), \
                patch.object(C.subprocess,'Popen',side_effect=launch):
            result=C.capture('/not-a-codec','0'*64,time.monotonic()+3,[Path(d)/'raw',Path(d)/'error'])
        self.assertEqual(result,(str(C.ROOT).encode()+b'\nabsent\n',b'',0,None))
        self.assertEqual(observed['cwd'],C.ROOT)
        self.assertEqual(observed['env'],dict(PATH='/usr/bin:/bin',LANG='C',LC_ALL='C',TZ='UTC',**C.SANITIZERS))

    def test_exact_pin_schema_and_shared_inputs(self):
        record=dict(path='/synthetic/pin',bytes=1,sha256='a'*64)
        self.assertEqual(C.R.pin_map([record,copy.deepcopy(record)]),{record['path']:record})
        for changed in (dict(record,bytes=True),dict(record,path='relative'),dict(record,sha256='A'*64)):
            with self.assertRaises(ValueError): C.R.pin_map([changed])
        with self.assertRaises(ValueError): C.R.pin_map([record,dict(record,bytes=2)])

    def test_archive_reader_reuse_is_data_only(self):
        self.assertIs(C.R.qualified_inputs,C.R.U.qualified_inputs)
        self.assertIsNot(C.R.build,C.R.U.build)
        source=(C.ROOT/'bench/Wh2K4CostBuildR0.py').read_text()
        for prohibited in ('U.build(', 'U.current(', 'U.replay(', 'U.run('):
            self.assertNotIn(prohibited,source)
        self.assertEqual(len(C.R.U.PRODUCERS),19)

    def test_new_builder_rejects_output_scope_before_work(self):
        with patch.object(C.R,'qualified_inputs',side_effect=AssertionError('unexpected provenance work')):
            for mode,path in [('wrong',Path('/tmp/native')),('native',C.ROOT/'native'),
                              ('native',Path('/var/tmp/native')),('native',Path('/tmp/other'))]:
                with self.assertRaises(ValueError): C.R.build(mode,path)

    def test_static_wrapper_and_shared_worker_unchanged(self):
        import hashlib
        wrapper=(C.ROOT/'bench/Wh2K4SerializedCostR0.cpp').read_text()
        self.assertIn('#define WH2_SMALL_COST_K 4\n',wrapper)
        self.assertIn('#define WH2_SMALL_COST_REPAIRS 8\n',wrapper)
        self.assertIn(C.PROTOCOL,wrapper)
        self.assertIn(str(C.OUTPUT/'CLAIM.json'),wrapper)
        baseline=C.command(['git','show','adb71c27257602ef372b16569d177ecae9ece776:bench/Wh2SmallLifecycleWorkerR0.h'])
        self.assertEqual(hashlib.sha256(baseline).digest(),
                         hashlib.sha256((C.ROOT/'bench/Wh2SmallLifecycleWorkerR0.h').read_bytes()).digest())


if __name__ == '__main__': unittest.main()
