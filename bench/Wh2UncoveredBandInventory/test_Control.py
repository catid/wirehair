"""Controller failure tests with fake artifacts and no native codec loading."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
import Control as C


class Test(unittest.TestCase):
    def test_spent_namespace_before_inputs_or_work(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            link=root/'dangling'; link.symlink_to(root/'absent')
            for path in (root,link):
                with mock.patch.object(C,'OUTPUT',path),mock.patch.object(C,'inputs') as inputs, \
                     mock.patch.object(C.Launch,'capture') as capture:
                    with self.assertRaisesRegex(ValueError,'namespace already spent'): C.run()
                    inputs.assert_not_called(); capture.assert_not_called()

    def test_bad_claim_and_interpreter_before_native_load(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); claim=root/'claim.json'
            claim.write_bytes(C.A.canonical({'interpreter':'/not-the-producing-python'}))
            with mock.patch.object(C,'OUTPUT',root),mock.patch.object(C.I,'run_inventory') as run, \
                 mock.patch.object(C,'current') as current:
                with self.assertRaisesRegex(ValueError,'actual inventory claim'): C.worker('0'*64)
                with self.assertRaisesRegex(ValueError,'claimed producing interpreter'):
                    C.worker(C.A.sha(claim.read_bytes()))
                run.assert_not_called(); current.assert_not_called()

    def test_replay_preserves_producing_interpreter(self):
        claim={'interpreter':'/recorded/producer','other':'unchanged'}
        with mock.patch.object(C,'inputs',return_value=claim) as inputs:
            C.current(claim)
            inputs.assert_called_once_with('/recorded/producer')
        with mock.patch.object(C,'inputs',return_value=dict(claim,other='changed')):
            with self.assertRaises(ValueError): C.current(claim)

    def test_capture_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            raw,err=(Path(directory)/p for p in ('raw','err'))
            raw.write_bytes(b'123'); err.write_bytes(b'')
            valid=dict(exit=0,failure=None,wall_seconds=1.0,stdout_bytes=3,stderr_bytes=0)
            C.check_capture(valid,raw,err)
            for key,value in (('exit',False),('exit',1),('failure','TIMEOUT'),
                              ('wall_seconds',float('nan')),('wall_seconds',float('inf')),
                              ('wall_seconds',0),('wall_seconds',151),('stdout_bytes',3.0),
                              ('stderr_bytes',False),('stdout_bytes',4),('extra',1)):
                with self.subTest(key=key,value=value):
                    with self.assertRaises(ValueError): C.check_capture(dict(valid,**{key:value}),raw,err)
            err.write_bytes(b'error')
            with self.assertRaises(ValueError): C.check_capture(valid,raw,err)

    def test_neutral_evidence_binds_source_and_producer(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); lib=root/'lib.so'; lib.write_bytes(b'fake')
            for name in C.RUNTIME: (root/name).write_bytes(name.encode())
            libraries={'native':(lib,C.A.sha(b'fake'))}
            with mock.patch.object(C,'HERE',root),mock.patch.object(C,'QUALIFICATION',root), \
                 mock.patch.object(C.I,'LIBRARIES',libraries), \
                 mock.patch.object(C.I,'verify',return_value={'parity_sha256':'same'}) as verify:
                paths=C.neutral_paths('native'); claim=C.neutral_inputs('native')
                C.A.publish(paths['claim'],claim)
                paths['raw'].write_bytes(b'fake'); paths['stderr'].write_bytes(b'')
                C.A.publish(paths['process'],dict(exit=0,failure=None,wall_seconds=1.0,
                                                stdout_bytes=4,stderr_bytes=0))
                result=C.neutral_evidence('native')
                self.assertEqual(result['artifacts']['claim'],C.A.pin(paths['claim']))
                verify.assert_called_once_with(paths['raw'],C.A.sha(paths['claim'].read_bytes()),'native',True)
                verify.reset_mock()
                (root/'Inventory.py').write_bytes(b'changed runtime')
                with self.assertRaisesRegex(ValueError,'current runtime sources'): C.neutral_evidence('native')
                verify.assert_not_called()

    def test_neutral_worker_rejects_changed_claim(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            with mock.patch.object(C,'QUALIFICATION',root), \
                 mock.patch.object(C,'neutral_inputs',return_value={'sources':'new'}), \
                 mock.patch.object(C.I,'run_inventory') as run:
                path=C.neutral_paths('native')['claim']; C.A.publish(path,{'sources':'old'})
                with self.assertRaisesRegex(ValueError,'unchanged neutral producer'):
                    C.neutral_worker('native',C.A.sha(path.read_bytes()))
                run.assert_not_called()

    def test_neutral_is_bounded_and_failed_outputs_are_not_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); claim={'interpreter':{'path':sys.executable}}
            def capture(command,out,err):
                self.assertEqual(command[3:5],['neutral-worker','native'])
                out.write_bytes(b'partial'); err.write_bytes(b'')
                return dict(exit=-9,failure='TIMEOUT',wall_seconds=150,stdout_bytes=7,stderr_bytes=0)
            with mock.patch.object(C,'QUALIFICATION',root),mock.patch.object(C,'codec_pins'), \
                 mock.patch.object(C,'neutral_inputs',return_value=claim), \
                 mock.patch.object(C.Launch,'capture',side_effect=capture) as captured, \
                 mock.patch.object(C.I,'run_inventory') as run, \
                 mock.patch.object(C,'neutral_evidence',side_effect=ValueError('failed capture')):
                with self.assertRaisesRegex(ValueError,'failed capture'): C.neutral('native')
                with self.assertRaisesRegex(ValueError,'fresh neutral outputs'): C.neutral('native')
                captured.assert_called_once(); run.assert_not_called()
                paths=C.neutral_paths('native')
                self.assertFalse(paths['proof'].exists())
                for key in ('claim','raw','stderr','process'):
                    self.assertEqual(paths[key].stat().st_mode & 0o777,0o400)

    def test_native_portable_parity_required(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); proofs={}
            with mock.patch.object(C,'QUALIFICATION',root):
                for mode in C.I.LIBRARIES:
                    paths=C.neutral_paths(mode)
                    claim=dict(sources={},interpreter=C.A.pin(sys.executable))
                    for name,path in paths.items():
                        if name not in ('proof','claim'): path.write_bytes(b'fake')
                    C.A.publish(paths['claim'],claim)
                    proofs[mode]=dict(result={'parity_sha256':mode})
                    C.A.publish(paths['proof'],proofs[mode])
                with mock.patch.object(C,'codec_pins',return_value={}), \
                     mock.patch.object(C,'neutral_evidence',side_effect=lambda mode:copy.deepcopy(proofs[mode])), \
                     mock.patch.object(C,'git') as git:
                    with self.assertRaisesRegex(ValueError,'all native/portable neutral records agree'): C.inputs()
                    git.assert_not_called()

    def test_scientific_failure_seals_spent_namespace(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)/'science'
            claim={'interpreter':sys.executable}
            with mock.patch.object(C,'OUTPUT',root),mock.patch.object(C,'inputs',return_value=claim), \
                 mock.patch.object(C.Launch,'capture',side_effect=OSError('synthetic launch failure')) as launch:
                with self.assertRaisesRegex(OSError,'synthetic launch failure'): C.run()
                failure=C.A.decode((root/'failed.json').read_bytes())
                self.assertEqual(failure['outcome'],'INVALID')
                self.assertTrue(failure['namespace_spent'])
                for path in root.iterdir(): self.assertEqual(path.stat().st_mode & 0o777,0o400)
                with self.assertRaisesRegex(ValueError,'namespace already spent'): C.run()
                launch.assert_called_once()

    def test_only_neutral_producer_can_launch_but_reader_can_differ(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); other=root/'reader'; other.write_bytes(b'another Python')
            producer=C.A.pin(sys.executable)
            proof=dict(result={'parity_sha256':'same'})
            with mock.patch.object(C,'QUALIFICATION',root),mock.patch.object(C,'OWN',()), \
                 mock.patch.object(C,'codec_pins',return_value={}),mock.patch.object(C,'head',return_value='fixed'), \
                 mock.patch.object(C,'neutral_evidence',return_value=proof):
                for name in ('python312.tests.log','python38.tests.log'): (root/name).write_bytes(b'fake')
                for mode in C.I.LIBRARIES:
                    paths=C.neutral_paths(mode)
                    C.A.publish(paths['claim'],dict(sources={},interpreter=producer))
                    C.A.publish(paths['proof'],proof)
                    for key in ('raw','stderr','process'): paths[key].write_bytes(b'fake')
                claim=C.inputs()
                with mock.patch.object(C.sys,'executable',str(other)):
                    with self.assertRaisesRegex(ValueError,'producer passed both neutral backends'): C.inputs()
                    C.current(claim)


if __name__=='__main__': unittest.main()
