#!/usr/bin/env python3
"""Neutral tests only: no scientific namespace or timed codec launch."""
import contextlib
import copy
import importlib.util
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

SPEC = importlib.util.spec_from_file_location('small_isolation_gate_tested',
    Path(__file__).with_name('Wh2SmallIsolationPreservedCostR0.py'))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
R = M.R


class Tests(unittest.TestCase):
    def test_explicit_configuration_preserves_original_defaults(self):
        default = R.configuration(None)
        self.assertNotEqual(M.PROTOCOL,R.PROTOCOL)
        self.assertNotEqual(M.OUTPUT,R.OUTPUT)
        self.assertEqual(default.libraries,R.N.LIBRARIES)
        self.assertNotEqual(default.libraries[1],M.LIBRARIES[1])
        self.assertEqual(M.SETTINGS.libraries[0],default.libraries[0])
        self.assertEqual(M.SETTINGS.sources[:3],R.NEW)
        self.assertEqual(len(R.CASES),20)
        self.assertEqual(len(list(R.roster())),51840)
        for cfg in (default,M.SETTINGS):
            result = R.combine([dict(outcome='PASS'),dict(outcome='REGRESSION')],cfg.protocol)
            self.assertEqual((result['protocol'],result['outcome']),(cfg.protocol,'REGRESSION'))
            self.assertFalse(result['WH1_speed_qualified'])

    def test_exact_native_compile_contract(self):
        relative = M.PREFIX+'codec/WirehairV2Profile.cpp.o'
        command = M.native_compile(relative)
        self.assertEqual(command[-2:],['-c',str(M.ROOT/'codec/WirehairV2Profile.cpp')])
        self.assertEqual(command[command.index('-o')+1],relative)
        self.assertEqual(command[command.index('-MF')+1],relative+'.d')
        self.assertIn('-O3',command)
        self.assertNotIn('-march=native',command)
        self.assertNotIn('-flto',command)

    def test_runtime_launcher_symlink_resolved_before_ldd(self):
        with tempfile.TemporaryDirectory(prefix='wh2-isolation-linkage-') as directory:
            folder=Path(directory); executable=folder/'actual'; launcher=folder/'launcher'
            library=folder/'library'; executable.write_bytes(b'neutral'); library.write_bytes(b'neutral')
            launcher.symlink_to(executable)
            with mock.patch.object(R,'command',return_value=('library => '+str(library)+' (0x1)\n').encode()) as cmd:
                self.assertEqual(R.runtime_dependencies(launcher),{library})
                cmd.assert_called_once_with(['/usr/bin/ldd',executable])

    def test_exact_candidate_elf_without_default_rebinding(self):
        before = R.metadata()
        actual = R.metadata(M.LIBRARIES)
        self.assertEqual(before[0],actual[0])
        self.assertNotEqual(before[1]['sha256'],actual[1]['sha256'])
        self.assertEqual(actual[1]['sha256'],M.LIBRARIES[1][1])
        self.assertEqual(R.metadata(),before)
        self.assertEqual([len(lib['exports']) for lib in actual],[53,53])
        self.assertEqual([len(lib['runtime_slots']) for lib in actual],[37,37])
        bad = (M.LIBRARIES[0],(M.LIBRARIES[1][0],'0'*64))
        with self.assertRaises(ValueError): R.metadata(bad)

    def test_dependency_roster_bounds_and_uniqueness(self):
        with tempfile.TemporaryDirectory(prefix='wh2-isolation-deps-') as directory:
            source = Path(directory)/'source'; source.write_bytes(b'neutral')
            relative = M.PREFIX+'codec/WirehairV2Profile.cpp.o'
            record = relative+': #deps 1, deps mtime 123 (VALID)\n    '+str(source)+'\n'
            self.assertEqual(M.dependency_roster(record),{relative:[str(source)]})
            self.assertEqual(M.dependency_roster('other.o: ignored\n\n'+record),{relative:[str(source)]})
            repeated=record.replace('#deps 1','#deps 2')+'    '+str(source)+'\n'
            self.assertEqual(M.dependency_roster(repeated),{relative:[str(source)]*2})
            for bad in (record.replace('(VALID)','(STALE)'),record.replace('#deps 1','#deps 2'),record+'\n'+record):
                with self.assertRaises(ValueError): M.dependency_roster(bad)

    def test_closed_candidate_source_object_chain_read_only(self):
        before = R.O.pin(M.BASE/'libwirehair.so.2.0.0')
        reports,inputs = M.provenance()
        self.assertEqual(len(reports),2)
        candidate = reports[1]
        self.assertEqual(candidate['source_head'],M.SOURCE_HEAD)
        self.assertEqual(candidate['changed_from_admission'],['WirehairV2Profile.cpp.o'])
        self.assertEqual(len(candidate['objects']),17)
        self.assertIn(M.ROOT/'codec/WirehairV2Profile.cpp',inputs)
        self.assertIn(M.BASE/'.ninja_deps',inputs)
        self.assertEqual(candidate['original'],before)
        self.assertEqual(R.O.pin(M.BASE/'libwirehair.so.2.0.0'),before)
        selected = [p for p in candidate['objects'] if Path(p['path']).name=='WirehairV2Profile.cpp.o']
        self.assertEqual(candidate['proof_object']['sha256'],selected[0]['sha256'])
        self.assertEqual(candidate['excluded_historical_diagnostics'][0]['claimed'],M.LOST_DIAGNOSTIC)
        self.assertNotEqual(candidate['excluded_historical_diagnostics'][0]['observed'],M.LOST_DIAGNOSTIC)

    def test_original_historical_gate_stays_strict(self):
        with self.assertRaises(ValueError): R.historical_claim(*R.HISTORICAL[1])
        path=R.N.LIBRARIES[1][0].parent/'libwirehair.a'
        forbidden={M.LOST_DIAGNOSTIC['path']:M.LOST_DIAGNOSTIC,str(path):R.O.pin(path)}
        with self.assertRaisesRegex(ValueError,'only explicitly recorded non-producing'):
            R.historical_claim(*R.HISTORICAL[1],diagnostic_exclusions=forbidden)

    def test_new_protocol_controller_replay_and_no_namespace_reuse(self):
        with tempfile.TemporaryDirectory(prefix='wh2-isolation-controller-') as directory:
            folder = Path(directory); receipt=folder/'receipt'; output=folder/'bundle'
            receipt.write_bytes(R.A.canonical(dict(executable='/not-a-codec')))
            cfg = M.SETTINGS._replace(output=output)
            launched=[]
            def capture(executable,claim,order,deadline,spools):
                self.assertEqual(executable,'/not-a-codec'); launched.append(order)
                R.A.publish(spools[0],b'raw\n'); R.A.publish(spools[1],b'')
                return b'raw\n',b'',0,None
            def current(frozen,settings):
                self.assertEqual(settings,cfg); return []
            def verify(raw,claim,order,meta,protocol):
                self.assertEqual(protocol,M.PROTOCOL); return dict(outcome='PASS')
            with mock.patch.object(R,'current',side_effect=current),mock.patch.object(R,'capture',side_effect=capture), \
                 mock.patch.object(R,'verify',side_effect=verify),contextlib.redirect_stdout(io.StringIO()):
                R.run(receipt,cfg)
                result = R.replay(cfg)
                self.assertEqual(result['protocol'],M.PROTOCOL)
                self.assertEqual(result['outcome'],'PASS')
                self.assertEqual(launched,[0,1])
                with self.assertRaises(FileExistsError): R.run(receipt,cfg)
                self.assertEqual(launched,[0,1])
                wrong = cfg._replace(protocol=R.PROTOCOL)
                with self.assertRaises(ValueError): R.replay(wrong)

    def test_current_requires_new_sources_and_recompiled_object(self):
        files={}; folder=Path('/synthetic/native'); exe=folder/'cost_worker'
        def add(path,data):
            files[str(path)]=data
            return dict(path=str(path),bytes=len(data),sha256=R.A.sha(data))
        def pin(path):
            data=files[str(path)]
            return dict(path=str(path),bytes=len(data),sha256=R.A.sha(data))
        inputs=[add(M.ROOT/name,b'source') for name in M.NEW]
        report=[dict(proof_name='proof-old.so',original=dict(sha256=R.A.sha(b'old'))),
                dict(proof_name='proof-new.so',original=dict(sha256=R.A.sha(b'new')),
                     proof_object=dict(name='proof-profile.o',sha256=R.A.sha(b'object')))]
        cfg = M.SETTINGS._replace(provenance=lambda: (report,set()))
        payloads={'cost_worker':b'worker','AdmissionLibraryBindings.h':b'bindings','library-metadata.json':b'[]',
                  'claim-binding.json':R.A.canonical(R.claim_binding(cfg)),'neutral-claim.json':b'neutral',
                  'library-provenance.json':R.A.canonical(report),'proof-old.so':b'old','proof-new.so':b'new',
                  'proof-profile.o':b'object','fixtures-old-new.json':b'{}','fixtures-new-old.json':b'{}',
                  'negative-cli.json':b'[]','link.map':b'link'}
        artifacts=[add(folder/name,data) for name,data in payloads.items()]
        manifest=dict(protocol=M.PROTOCOL,mode='native',scientific_launch=False,library_source_provenance_closed=True,
                      sanitized_library_code=False,environment={k:None for k in R.ENV_KEYS+('ASAN_OPTIONS','UBSAN_OPTIONS')},
                      inputs=inputs,artifacts=artifacts)
        def frozen():
            mp=add(folder/'manifest.json',R.A.canonical(manifest))
            return dict(protocol=M.PROTOCOL,head='a'*40,executable=str(exe),
                        environment={k:None for k in R.ENV_KEYS},
                        pins=sorted(manifest['inputs']+manifest['artifacts']+[mp],key=lambda p:p['path']))
        with mock.patch.dict(os.environ,{},clear=True), \
             mock.patch.object(R,'command',side_effect=lambda args:b'a'*40+b'\n' if args[1]=='rev-parse' else b'source'), \
             mock.patch.object(R.O,'pin',side_effect=pin), \
             mock.patch.object(R.A,'read_regular',side_effect=lambda p,cap:files[str(p)]), \
             mock.patch.object(R,'metadata',return_value=[]), \
             mock.patch.object(R,'bindings_header',return_value=b'bindings'),mock.patch.object(R,'verify_header'):
            R.current(frozen(),cfg)
            saved=copy.deepcopy(manifest)
            for missing in (str(folder/'proof-profile.o'),str(M.ROOT/M.NEW[-1]),str(M.ROOT/M.NEW[-2])):
                manifest=copy.deepcopy(saved)
                for key in ('inputs','artifacts'):
                    manifest[key]=[p for p in manifest[key] if p['path']!=missing]
                with self.assertRaises(ValueError): R.current(frozen(),cfg)
            manifest=copy.deepcopy(saved)
            files[str(folder/'proof-profile.o')]=b'wrong'
            with self.assertRaises(ValueError): R.current(frozen(),cfg)


if __name__=='__main__': unittest.main()
