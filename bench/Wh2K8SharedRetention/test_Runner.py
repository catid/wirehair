"""Synthetic controller/receipt tests; never execute a codec or real cohort."""
import contextlib
import copy
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

import Runner as C

K, R, A = C.K, C.K.R, C.K.A


class RunnerTest(unittest.TestCase):
    def test_configuration_keeps_exact_observer_and_adds_mandatory_sources(self):
        old, cfg = K.settings(C.ADAPTER), C.settings()
        self.assertEqual(cfg._replace(sources=old.sources, provenance=old.provenance), old)
        self.assertEqual(cfg.sources, old.sources+C.NEW)
        self.assertIs(cfg.provenance, C.provenance)
        with patch.object(K, 'provenance', return_value=(['new-proof'], {Path('/library')})) as source, \
             patch.object(C, 'observer_inputs', return_value={Path('/neutral')}):
            self.assertEqual(C.provenance(Path('/new-build')), (['new-proof'], {Path('/library'),Path('/neutral')}))
            source.assert_called_once_with(C.ADAPTER, Path('/new-build'))

    def test_neutral_manifest_closure_rejects_drift(self):
        files = {}
        def add(path, raw):
            files[str(path)] = raw
            return dict(path=str(path), bytes=len(raw), sha256=A.sha(raw))
        shared = add('/synthetic/shared', b'shared')
        expected, paths = [], {Path(shared['path'])}
        for mode,_ in C.MANIFESTS:
            source = add('/synthetic/'+mode, mode.encode())
            manifest = dict(protocol=K.PROTOCOL, mode=mode, scientific_launch=False,
                            library_source_provenance_closed=True, sanitized_library_code=False,
                            inputs=[shared], artifacts=[source])
            path = C.QUALIFIED/mode/'manifest.json'
            pin = add(path, A.canonical(manifest)); expected.append((mode,pin['sha256']))
            paths.update((path, Path(source['path'])))
        def pin(path):
            raw=files[str(path)]; return dict(path=str(path), bytes=len(raw), sha256=A.sha(raw))
        with patch.object(C, 'MANIFESTS', tuple(expected)), patch.object(K.B, 'pin', side_effect=pin), \
             patch.object(A, 'read_regular', side_effect=lambda path,cap:files[str(path)]):
            self.assertEqual(C.observer_inputs(), paths)
            files['/synthetic/shared'] = b'changed'
            with self.assertRaisesRegex(ValueError, 'unchanged neutral prerequisite'): C.observer_inputs()
            files['/synthetic/shared'] = b'shared'
            files[str(C.QUALIFIED/'asan-driver/manifest.json')] = b'{}'
            with self.assertRaisesRegex(ValueError, 'accepted neutral observer manifest'): C.observer_inputs()

    def test_main_always_supplies_new_configuration(self):
        cfg = C.settings()
        with patch.object(C, 'enter_clean_environment'), patch.object(C, 'settings', return_value=cfg), \
             patch.object(K.S, 'build') as build, patch.object(R, 'receipt', return_value={'receipt':'fake'}) as receipt, \
             patch.object(R, 'run') as run, patch.object(R, 'replay', return_value={'outcome':'PASS'}) as replay, \
             patch.object(A, 'publish') as publish, contextlib.redirect_stdout(io.StringIO()):
            C.main(['build','native','/new/native'])
            build.assert_called_once_with('native',Path('/new/native'),cfg,R)
            C.main(['receipt','/new/native','/new/receipt'])
            receipt.assert_called_once_with(Path('/new/native'),cfg)
            publish.assert_called_once_with(Path('/new/receipt'),A.canonical({'receipt':'fake'}))
            C.main(['run','/new/receipt']); run.assert_called_once_with(Path('/new/receipt'),cfg)
            C.main(['replay']); replay.assert_called_once_with(cfg)

    def test_environment_reexec_removes_unlisted_policy(self):
        environment = K.B.process_environment()
        with patch.dict(os.environ,dict(environment,LD_BIND_NOT='1',MALLOC_ARENA_MAX='1'),clear=True), \
             patch.object(C.os,'execve',side_effect=RuntimeError('captured exec')) as execute:
            with self.assertRaisesRegex(RuntimeError,'captured exec'): C.enter_clean_environment()
            self.assertEqual(execute.call_args[0][0],C.sys.executable)
            self.assertEqual(execute.call_args[0][2],environment)
        with patch.dict(os.environ,environment,clear=True), patch.object(C.os,'execve') as execute:
            C.enter_clean_environment(); execute.assert_not_called()

    def test_current_receipt_cannot_drop_pins_rebind_or_use_old_protocol(self):
        files = {}
        def add(path,raw):
            files[str(path)] = raw
            return dict(path=str(path),bytes=len(raw),sha256=A.sha(raw))
        cfg = C.settings()._replace(header_checker=Mock())
        folder=Path('/synthetic/native'); exe=folder/'cost_worker'
        inputs=[add(path,b'source') for path in sorted({K.ROOT/name for name in cfg.sources})]
        inputs.append(add('/synthetic/prerequisite',b'neutral'))
        provenance=[dict(proof_name='proof-old.so',original=dict(sha256=A.sha(b'old'))),
                    dict(proof_name='proof-new.so',original=dict(sha256=A.sha(b'new')))]
        cfg=cfg._replace(provenance=Mock(return_value=(provenance,{Path('/synthetic/prerequisite')})))
        payloads={'cost_worker':b'worker','AdmissionLibraryBindings.h':b'bindings','library-metadata.json':b'[]',
                  'library-provenance.json':A.canonical(provenance),'proof-old.so':b'old','proof-new.so':b'new',
                  'claim-binding.json':A.canonical(R.claim_binding(cfg)),'neutral-claim.json':b'neutral',
                  'fixtures-old-new.json':b'{}','fixtures-new-old.json':b'{}','negative-cli.json':b'[]','link.map':b'link'}
        artifacts=[add(folder/name,raw) for name,raw in payloads.items()]
        manifest=dict(protocol=cfg.protocol,mode='native',scientific_launch=False,library_source_provenance_closed=True,
                      sanitized_library_code=False,environment={key:None for key in R.ENV_KEYS+('ASAN_OPTIONS','UBSAN_OPTIONS')},
                      inputs=inputs,artifacts=artifacts)
        manifest_pin=add(folder/'manifest.json',A.canonical(manifest))
        frozen=dict(protocol=cfg.protocol,head='a'*40,executable=str(exe),environment={key:None for key in R.ENV_KEYS},
                    pins=sorted(inputs+artifacts+[manifest_pin],key=lambda pin:pin['path']))
        def pin(path):
            raw=files[str(path)]; return dict(path=str(path),bytes=len(raw),sha256=A.sha(raw))
        with patch.dict(os.environ,{},clear=True), patch.object(R,'command',side_effect=lambda args:b'a'*40+b'\n' if args[1]=='rev-parse' else b'source'), \
             patch.object(R.O,'pin',side_effect=pin), patch.object(R.A,'read_regular',side_effect=lambda path,cap:files[str(path)]), \
             patch.object(R,'metadata',return_value=[]) as metadata, patch.object(R,'bindings_header',return_value=b'bindings'):
            R.current(frozen,cfg); metadata.assert_called_once_with(K.LIBRARIES)
            for i in range(len(frozen['pins'])):
                bad=copy.deepcopy(frozen); del bad['pins'][i]
                with self.assertRaises(ValueError): R.current(bad,cfg)
            bad=copy.deepcopy(frozen); bad['pins'].append(bad['pins'][0])
            with self.assertRaises(ValueError): R.current(bad,cfg)
            for field,value in (('protocol',R.PROTOCOL),('head','b'*40),
                                ('executable','/synthetic/asan-driver/cost_worker')):
                bad=copy.deepcopy(frozen); bad[field]=value
                with self.assertRaises(ValueError): R.current(bad,cfg)

    def test_invalid_controller_keeps_both_orders_and_cannot_retry(self):
        for failure_kind in ('capture','validation'):
            with tempfile.TemporaryDirectory(prefix='wh2-k8-controller-test.') as directory:
                root=Path(directory); cfg=C.settings()._replace(output=root/'bundle')
                receipt=root/'receipt'; A.publish(receipt,A.canonical(dict(executable='/not-a-codec')))
                launched=[]
                def capture(executable,claim,order,deadline,spools):
                    self.assertEqual(executable,'/not-a-codec'); launched.append(order)
                    A.publish(spools[0],b'prefix\n'); A.publish(spools[1],b'')
                    return b'prefix\n',b'',0,('synthetic capture error' if order==0 and failure_kind=='capture' else None)
                verify_results=[dict(outcome='PASS')] if failure_kind=='capture' else [ValueError('bad first order'),dict(outcome='PASS')]
                with patch.object(R,'current',return_value=[]), patch.object(R,'capture',side_effect=capture), \
                     patch.object(R,'verify',side_effect=verify_results), contextlib.redirect_stdout(io.StringIO()):
                    R.run(receipt,cfg)
                    self.assertEqual(launched,[0,1])
                    complete=A.decode((cfg.output/'COMPLETE.json').read_bytes())
                    self.assertEqual(complete['outcome'],'INVALID')
                    self.assertEqual(len(A.decode((cfg.output/'processes.json').read_bytes())),2)
                    with self.assertRaises(FileExistsError): R.run(receipt,cfg)
                    self.assertEqual(launched,[0,1])
                    self.assertTrue(all(Path(pin['path']).stat().st_mode&0o777==0o400 for pin in complete['files']))

    def test_preflight_failure_does_not_claim_namespace(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k8-preflight-test.') as directory:
            root=Path(directory); cfg=C.settings()._replace(output=root/'bundle')
            receipt=root/'receipt'; A.publish(receipt,A.canonical(dict(executable='/not-a-codec')))
            with patch.object(R,'current',side_effect=ValueError('bad prerequisite')), patch.object(R,'capture') as capture:
                with self.assertRaisesRegex(ValueError,'bad prerequisite'): R.run(receipt,cfg)
                capture.assert_not_called(); self.assertFalse(cfg.output.exists())

    def test_controller_and_replay_use_exact_new_combiner(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k8-replay-test.') as directory:
            root=Path(directory); cfg=C.settings()._replace(output=root/'bundle')
            receipt=root/'receipt'; A.publish(receipt,A.canonical(dict(executable='/not-a-codec')))
            def capture(executable,claim,order,deadline,spools):
                A.publish(spools[0],b'raw\n'); A.publish(spools[1],b'')
                return b'raw\n',b'',0,None
            with patch.object(R,'current',return_value=[]), patch.object(R,'capture',side_effect=capture), \
                 patch.object(R,'verify',return_value=dict(outcome='PASS')) as verify, contextlib.redirect_stdout(io.StringIO()):
                R.run(receipt,cfg); result=R.replay(cfg)
                self.assertTrue(result['current_path_retention_qualified'])
                self.assertFalse(result['pre_admission_restoration_qualified'])
                self.assertFalse(result['historical_K3_workload_retention_qualified'])
                self.assertNotIn('certified_four_cell_improvements',result)
                for call in verify.call_args_list:
                    self.assertEqual(call[0][4:6],(K.PROTOCOL,K.CASES))
                    self.assertIs(call[0][6],K.verify_header)
                A.publish(cfg.output/'unexpected',b'extra')
                with self.assertRaises(ValueError): R.replay(cfg)


if __name__ == '__main__': unittest.main()
