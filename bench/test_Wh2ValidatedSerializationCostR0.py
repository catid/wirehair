"""Synthetic and source-only checks; never launch a codec or timing cohort."""
import copy
import importlib.util
from pathlib import Path
import types
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('serialization_tested',
    Path(__file__).with_name('Wh2ValidatedSerializationCostR0.py'))
C = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(C)
B = C.B

def synthetic(ratio=1):
    records = []
    for c in C.R.roster(C.CASES):
        duration = round(1000000*(ratio if c[5]==2 and C.R.SIDES[c[6]]^c[2] else 1))
        records.append(dict(coordinate=c,observation=dict(clocks=[0,0,0,duration,0,0])))
    return C.R.statistics(records,C.CASES)

class Tests(unittest.TestCase):
    def test_historical_defaults_unchanged(self):
        cfg = C.R.configuration(None)
        self.assertEqual(cfg.cases,C.R.CASES)
        self.assertIsNone(cfg.header_checker); self.assertIsNone(cfg.result_combiner)
        self.assertEqual(C.R.callback_count(),51840)
        self.assertEqual(C.D.SETTINGS.cases,C.R.CASES)

    def test_complete_roster_and_phases(self):
        self.assertEqual(C.CASES[:20],C.R.CASES)
        self.assertEqual(len(C.CASES),38)
        self.assertEqual(len(set(C.CASES)),38)
        self.assertEqual(C.R.callback_count(C.CASES),98496)
        seen,phases = {},{}
        for index,c in enumerate(C.R.roster(C.CASES)):
            i,r,o,w,m,pair,p,a,q = c
            self.assertEqual(i,index)
            self.assertEqual(a,C.R.PAIRS[pair][C.R.SIDES[p]^o])
            key = w,m,pair,o
            if p==0: seen.setdefault(key,[]).append(r)
            if p>=2 and p%2==0: phases.setdefault(key,[]).append((q*96//1000000)//2)
        self.assertEqual(len(seen),456)
        self.assertEqual(sum(k[2]<2 for k in seen),304)
        for key in seen:
            self.assertEqual(sorted(seen[key]),list(range(12)))
            self.assertEqual(sorted(phases[key]),sorted(list(range(48))*2))

    def test_small_rows_against_separate_oracles(self):
        k3 = C.R.O.selected_rows()
        self.assertEqual(C.rows(3),k3[:11]+k3[12:])
        for k in (5,8):
            oracle = C.R.sibling('serialization_oracle_'+str(k),'Wh2K%dPublicRecoveryR0.py'%k)
            self.assertEqual(C.rows(k),tuple(oracle.coefficient(i) for i in C.packet_ids(k)))
        for k in (3,5,8):
            self.assertEqual(C.rows(k)[:k],tuple(tuple(int(i==j) for j in range(k)) for i in range(k)))
            for b in (2,64,1280):
                arm = C.small_fixture(k,b)
                self.assertEqual(len(bytes.fromhex(arm['packets'])),(k+14)*b)
                self.assertGreaterEqual(arm['steps'],k)
                self.assertEqual(bytes.fromhex(arm['packets'])[:k*b],bytes((37*i+i//11)%256 for i in range(k*b)))
                self.assertEqual(len(bytes.fromhex(arm['profile'])),32)
        self.assertIn(12,C.packet_ids(5))
        self.assertEqual(len(C.packet_ids(8)),22)
        with self.assertRaises(ValueError): C.rows(6)

    def test_rank_uses_actual_combined_order(self):
        self.assertEqual(C.first_success([(1,0),(2,0),(0,1)],2),3)
        self.assertEqual(C.first_success([(0,0),(0,1),(1,0)],2),3)
        with self.assertRaises(ValueError): C.first_success([(1,0),(2,0)],2)
        for k in (3,5,8):
            self.assertEqual(C.small_fixture(k,2)['steps'],C.first_success(C.rows(k)[k:]+C.rows(k)[:k],k))

    def test_small_fixture_tampering(self):
        fixtures = [None]*20
        for c in C.CASES[20:]:
            _,k,b,_ = c
            fixtures.append(dict(case=list(c),batch=128,
                source=bytes((37*i+i//11)%256 for i in range(k*b)).hex(),arms=[C.small_fixture(k,b)]*2))
        header = dict(fixtures=fixtures)
        with patch.object(C.R,'verify_header') as old:
            C.verify_header(header,0,'claim',[])
            self.assertEqual(len(old.call_args[0][0]['fixtures']),20)
            for slot in range(20,38):
                for key,value in (('case',[4,8,1280,2]),('batch',4),('source',''),('arms',[])):
                    bad = copy.deepcopy(header); bad['fixtures'][slot][key] = value
                    if bad==header: continue
                    with self.assertRaises(ValueError): C.verify_header(bad,0,'claim',[])
                for field in ('profile','packets','steps'):
                    bad = copy.deepcopy(header); bad['fixtures'][slot]['arms'][1][field] = 0
                    with self.assertRaises(ValueError): C.verify_header(bad,0,'claim',[])
            for changed in (fixtures[:-1],fixtures+[fixtures[-1]]):
                with self.assertRaises(ValueError): C.verify_header(dict(fixtures=changed),0,'claim',[])

    def test_retention_requires_benefit_and_both_orders(self):
        same,better = synthetic(),synthetic(.99)
        self.assertEqual(len(same['statistics']),456)
        self.assertEqual(C.combine([same,same])['outcome'],'PASS')
        self.assertFalse(C.combine([same,same])['candidate_retained'])
        self.assertTrue(C.combine([better,better])['candidate_retained'])
        self.assertFalse(C.combine([better,same])['candidate_retained'])
        self.assertFalse(C.combine([same,better])['candidate_retained'])
        for outcome in ('CONTROL_FAIL','REGRESSION','INCONCLUSIVE'):
            failed = copy.deepcopy(better); failed['outcome'] = outcome
            self.assertEqual(C.combine([better,failed])['outcome'],outcome)
            self.assertFalse(C.combine([better,failed])['candidate_retained'])
        for key in ('WH1_speed_qualified','all_K_claimed','static_speed_qualified',
                    'recovery_rate_claimed','production_promotion_claimed','pre_admission_restoration_qualified'):
            self.assertFalse(C.combine([better,better])[key])

    def test_no_five_percent_floor_and_no_resolved_regression(self):
        slight = synthetic(.999)
        self.assertTrue(C.combine([slight,slight])['candidate_retained'])
        slower = synthetic(1.001)
        self.assertEqual(slower['outcome'],'REGRESSION')
        self.assertFalse(C.combine([slower,slower])['candidate_retained'])

    def test_full_raw_checks_new_k8_cases_and_footer(self):
        zero = [0]*4
        previous = dict(clocks=[1,0,2,3,0,4],before=zero,after=zero)
        header = dict(prelude=previous,fixtures=[dict(arms=[dict(steps=c[1])]*2) for c in C.CASES])
        raw = [header]; addresses = {n:[123]*n+[0]*(128-n) for n in (4,128)}
        for c in C.R.roster(C.CASES):
            ready = previous['clocks'][5]+1; cpu = previous['clocks'][4]
            target = ready+c[-1]; start = target+3; case = C.CASES[c[3]]; cycles = C.R.batch(case)
            k,metric = case[1],c[4]
            observed = dict(clocks=[start,cpu+2,start+1,start+1000001,cpu+1000002,start+1000002],before=zero,after=zero)
            raw.append(dict(type='record',coordinate=c,ready=ready,target=target,
                wait=[ready+1,cpu,target+2,cpu+1],observation=observed,
                counts=[0 if metric else cycles,0 if metric else cycles*(k+14),cycles if metric else 0,
                        cycles*k if metric else 0,cycles if metric else 0,cycles],
                addresses=addresses[cycles],address_count=cycles,complete=True,checked=True))
            previous = observed
        count = C.R.callback_count(C.CASES)
        raw.append(dict(type='footer',complete=True,records=count,work_ns=count*1000000))
        with patch.object(C,'verify_header') as checker:
            self.assertEqual(C.R.verify_rows(raw,'claim',1,[],C.PROTOCOL,C.CASES,checker)['outcome'],'PASS')
            checker.assert_called_once_with(header,1,'claim',[],C.PROTOCOL)
            index = next(i for i,r in enumerate(raw[1:-1],1) if r['coordinate'][3]==37)
            saved = raw[index]; raw[index] = dict(saved,counts=[0]*6)
            with self.assertRaisesRegex(ValueError,'every attempted API call'):
                C.R.verify_rows(raw,'claim',1,[],C.PROTOCOL,C.CASES,checker)
            raw[index] = saved; raw[-1]['records'] = 82944
            with self.assertRaisesRegex(ValueError,'terminal complete footer'):
                C.R.verify_rows(raw,'claim',1,[],C.PROTOCOL,C.CASES,checker)
            with self.assertRaisesRegex(ValueError,'whole raw cohort'):
                C.R.verify_rows(raw,'claim',1,[])

    def test_generated_worker_isolation_and_limits(self):
        common,worker = (b.decode() for b in B.worker_sources())
        original = (C.ROOT/'bench/Wh2AdmissionRegressionCostR0.cpp').read_text()
        def work(text): return text[text.index('NOINLINE void RunWork('):text.index('void Prepare(')]
        self.assertEqual(work(common),work(original))
        self.assertTrue(worker.startswith('#define WH2_ADMISSION_PUBLIC_SMALL 1\n'))
        self.assertIn('#include "Wh2ValidatedSerializationCommon.h"',worker)
        self.assertIn('case_count=38,max_batch=128,callbacks=98496',common)
        self.assertIn('max_batch*(22*1280+128)',common)
        self.assertEqual(common.count('WIREHAIR_V2_PROFILE_SMALL_K8_2026_09'),2)
        for text in (common,worker):
            self.assertIn('cpu={210,210}',text)
            self.assertNotIn('UINT64_C(180000000000)',text)
            self.assertIn('UINT64_C(150000000000)',text)
        with self.assertRaises(ValueError): B.exact_replace('aaa','a','b')
        with self.assertRaises(ValueError): B.exact_replace('aaa','x','b')

    def test_shared_and_observer_link_inputs(self):
        self.assertTrue({'crtbeginS.o','crtendS.o','libdl.a'}<=set(B.LINK_INPUTS))
        command = B.link(Path('/tmp/new.so'),[Path('/tmp/a.o'),Path('/tmp/b.o')],Path('/tmp/new.map'))
        self.assertEqual(command[-4:],['/tmp/a.o','/tmp/b.o','-lm','-Wl,-Map,/tmp/new.map'])
        self.assertIn('-shared',command)
        self.assertIn('-Wl,-soname,libwirehair.so.2',command)
        files = {C.ROOT/'codec/WirehairV2Profile.cpp',C.ROOT/'gf256.cpp'}
        raw = ''.join('LOAD '+str(p)+'\n' for p in sorted(files)).encode()
        with patch.object(B.A,'read_regular',return_value=raw):
            self.assertEqual(B.verify_link_inputs(Path('/tmp/mock.map'),files),files)
            with self.assertRaisesRegex(ValueError,'frozen before link'):
                B.verify_link_inputs(Path('/tmp/mock.map'),set())
        for raw in (b'',b'LOAD relative.o\n'):
            with patch.object(B.A,'read_regular',return_value=raw):
                with self.assertRaises(ValueError): B.verify_link_inputs(Path('/tmp/mock.map'),files)

    def test_unregistered_import_closure(self):
        parent = types.ModuleType('fake_parent'); child = types.ModuleType('fake_child')
        parent.__file__ = str(C.ROOT/'bench/Wh2ValidatedSerializationCostBuildR0.py')
        child.__file__ = str(C.ROOT/'bench/Wh2AdmissionRegressionNeutral.py')
        parent.child = child; child.parent = parent
        with patch.object(B,'P',parent):
            found = B.imported_files()
        self.assertIn(Path(child.__file__),found)
        self.assertIn(Path(parent.__file__),found)

    def test_preparation_pin_is_mandatory(self):
        with self.assertRaises(ValueError): B.prepared_inputs('')
        with patch.object(B,'pin',return_value=dict(sha256='a'*64)):
            with self.assertRaisesRegex(ValueError,'exact independently reviewed preparation'):
                B.prepared_inputs('b'*64)
        with patch.object(C,'CANDIDATE_SHA',''):
            with self.assertRaises(ValueError): C.settings()

    def test_mutable_exceptions_cannot_include_production(self):
        self.assertEqual(set(B.MUTABLE_HARNESS),{
            'bench/Wh2ValidatedSerializationCostR0.py','bench/Wh2ValidatedSerializationCostBuildR0.py',
            'bench/test_Wh2ValidatedSerializationCostR0.py','bench/Wh2ValidatedSerializationCostR0.md',
            'bench/Wh2ValidatedSerializationRuntimeR0.py'})
        self.assertTrue(set(B.MUTABLE_HARNESS)<=set(C.SOURCES))
        self.assertNotIn('codec/WirehairV2Profile.cpp',B.MUTABLE_HARNESS)

    def test_controller_reexec_removes_unlisted_loader_policy(self):
        environment = B.P.process_environment()
        with patch.dict(C.os.environ,dict(environment,LD_BIND_NOT='1',MALLOC_ARENA_MAX='1'),clear=True):
            with patch.object(C.os,'execve',side_effect=RuntimeError('captured exec')) as execute:
                with self.assertRaisesRegex(RuntimeError,'captured exec'): C.enter_clean_environment()
                self.assertEqual(execute.call_args[0][2],environment)
                self.assertEqual(execute.call_args[0][0],C.sys.executable)
        with patch.dict(C.os.environ,environment,clear=True):
            with patch.object(C.os,'execve') as execute:
                C.enter_clean_environment(); execute.assert_not_called()

    def test_source_freezes_worker_before_qualification(self):
        text = Path(B.__file__).read_text()
        linked = text.index("run(args+list(map(str,objects))")
        freeze = text.index('dependencies.add(executable); P.freeze_inputs(dependencies,frozen)',linked)
        first_use = text.index("run(['/usr/bin/nm','-g',executable])",linked)
        self.assertLess(linked,freeze); self.assertLess(freeze,first_use)
        self.assertIn('P.freeze_inputs(dependencies,frozen)\n    cfg.qualify(',text)

    def test_preparation_order_and_replacement_mutations(self):
        objects = [str(B.BASE/'CMakeFiles/wirehair_objects.dir'/(name+'.o')) for name in B.P.PRODUCERS]
        replacement = [str(B.PREPARED/'proof-profile.o') if Path(p).name=='WirehairV2Profile.cpp.o' else p for p in objects]
        old = dict(path=str(B.DSO),bytes=686992,sha256=B.DSO_SHA)
        candidate = dict(path=str(B.PREPARED/'libwirehair.so.2.0.0'),bytes=1,sha256='b'*64)
        report = dict(original=old,candidate=candidate,original_objects=objects,candidate_objects=replacement,
            commands=[],inputs=[],artifacts=[],scientific_launch=False,producing_source_closure=True)
        raw = B.A.canonical(report)
        def pin(path):
            return candidate if Path(path).name=='libwirehair.so.2.0.0' else dict(sha256='a'*64)
        with patch.object(B,'pin',side_effect=pin), patch.object(B.A,'read_regular',return_value=raw), \
             patch.object(B,'worker_sources',return_value=()), patch.object(B,'overlay_inputs',return_value=set()), \
             patch.object(B.T,'runtime_source',return_value=raw):
            self.assertEqual(B.prepared_inputs('a'*64)[0],report)
            for field,value in (('original_objects',objects[::-1]),('candidate_objects',objects),
                                ('candidate_objects',replacement[:-1]),('scientific_launch',True),
                                ('producing_source_closure',False)):
                bad = copy.deepcopy(report); bad[field] = value
                with patch.object(B.A,'read_regular',return_value=B.A.canonical(bad)):
                    with self.assertRaises(ValueError): B.prepared_inputs('a'*64)

    def test_explicit_slot_rosters(self):
        for index, names in enumerate((B.T.SLOTS,B.T.SLOTS[:-1])):
            good = types.SimpleNamespace(slots=[(n,8*(i+1)) for i,n in enumerate(names)],
                                         exports=set(B.T.SLOTS))
            C.R.validate_slot_roster(good,index)
            for change in ('wrong-name','duplicate-name','missing','extra','duplicate-offset','missing-export'):
                bad = copy.deepcopy(good)
                if change=='wrong-name': bad.slots[0] = ('wirehair_wrong',bad.slots[0][1])
                elif change=='duplicate-name': bad.slots[0] = (bad.slots[1][0],bad.slots[0][1])
                elif change=='missing': bad.slots.pop()
                elif change=='extra': bad.slots.append(('wirehair_extra',4096))
                elif change=='duplicate-offset': bad.slots[0] = (bad.slots[0][0],bad.slots[1][1])
                else: bad.exports.remove('wirehair_v2_profile_serialize')
                with self.assertRaises(ValueError): C.R.validate_slot_roster(bad,index)
            with self.assertRaises(ValueError): C.R.validate_slot_roster(good,1-index)
        for index in (-1,2,True):
            with self.assertRaises(ValueError): C.R.validate_slot_roster(good,index)

    def test_metadata_identity_and_private_globals(self):
        for f in (C.R.current,C.R.receipt,C.R.run,C.R.replay):
            self.assertIs(f.__globals__,vars(C.R))
            self.assertIs(f.__globals__['metadata'],C.R.metadata)
            self.assertIs(f.__globals__['bindings_header'],C.R.bindings_header)
        self.assertIsNot(C.R,C.D.R)
        self.assertIn("A.exact(len(elf.slots),6,'exact internal public GOT roster')",
                      Path(C.D.R.__file__).read_text())
        right = ((B.DSO,B.DSO_SHA),(B.PREPARED/'libwirehair.so.2.0.0','a'*64))
        for libs in (None,(),right[:1],right[::-1],right+right,((right[0][0],'0'*64),right[1])):
            with self.assertRaises(ValueError): C.R.metadata(libs)
        with self.assertRaises(ValueError): B.T.runtime(Path('/tmp/not-the-derived-runtime.py'))
        self.assertEqual(B.T.PREPARED,B.PREPARED)
        self.assertEqual(B.T.BASELINE,B.DSO)

    def test_counted_cpp_and_header(self):
        common,_ = (v.decode() for v in B.worker_sources())
        self.assertNotIn('for(const auto& s:spec.slots)',common)
        self.assertEqual(common.count('i<spec.slot_count'),2)
        self.assertEqual(common.count('spec.slot_count<=spec.slots.size()'),2)
        self.assertIn('slot_count==(index==0?6u:5u)',common)
        self.assertIn('slot_count==(a==0?6u:5u)',common)
        meta = []
        for names in (B.T.SLOTS,B.T.SLOTS[:-1]):
            meta.append(dict(path='mock',sha256='a'*64,exports=[],
                slots=[dict(name=n,offset=i*8) for i,n in enumerate(names)],
                runtime_slots=[],context=100,context_bytes=200,getter=300))
        header = C.R.bindings_header(meta).decode()
        self.assertIn('}},6,{{}},100,200,300}',header)
        self.assertIn('}},5,{{}},100,200,300}',header)
        self.assertEqual(header.count('wirehair_v2_profile_serialize'),1)

    def test_actual_readonly_metadata_when_prepared(self):
        if not C.CANDIDATE_SHA: self.skipTest('library preparation not yet pinned')
        meta = C.R.metadata(C.settings().libraries)
        self.assertEqual([len(m['slots']) for m in meta],[6,5])
        self.assertEqual([len(m['exports']) for m in meta],[53,53])
        self.assertEqual([len(m['runtime_slots']) for m in meta],[37,37])

if __name__=='__main__': unittest.main()
