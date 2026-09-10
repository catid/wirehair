"""Neutral tests only: unrelated matrices and synthetic K8 outcomes."""
import hashlib
import importlib.util
import itertools
from pathlib import Path
import tempfile
import types
import unittest
from unittest import mock


def sibling(name,filename):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(filename))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


M=sibling('k8_tm_tested','Wh2K8ThueMorseR0.py');R=M.R
OLD=sibling('k8_controller_tests','test_Wh2NoncommutingRadixRunR0.py');OLD.M=R.C
NEUTRAL=(5,7,11,13,17,19,23,29)


class Tests(unittest.TestCase):
    def test_field_and_rank(self):
        M.F.init_field();M.G.reference_rank([])
        self.assertEqual(M.G.REFERENCE[0],M.F.MUL)
        self.assertEqual(M.G.REFERENCE[1],M.F.INV)
        for i in range(256):
            raw=hashlib.sha512(('k8-neutral/'+str(i)).encode()).digest()
            rows=[raw[j:j+8] for j in range(0,64,8)]
            self.assertEqual(M.checked_rank(rows),M.G.reference_rank(rows))
        for rows in ([],[(0,)*8],M.K.identity(8),[(1,)+(0,)*7]*12):
            self.assertEqual(M.checked_rank(rows),M.G.reference_rank(rows))
        with mock.patch.object(M.F,'matrix_rank',return_value=0):
            with self.assertRaisesRegex(ValueError,'rank disagreement'):
                M.checked_rank(M.K.identity(8))

    def test_unrelated_dimension8_mapper(self):
        pair=(M.K.companion(NEUTRAL),M.K.companion((9,)+NEUTRAL[1:]))
        budget=M.F.Budget();mapper=M.K.Mapper(pair,budget)
        evidence=M.verify_mapper(pair,mapper,budget)
        self.assertEqual(evidence['lookup_bytes'],65536)
        self.assertEqual(evidence['lookup_sha256'],hashlib.sha256(mapper.payload).hexdigest())
        for bit in range(2,32):
            for offset in (-1,0,1):
                i=(1<<bit)+offset
                self.assertEqual(mapper.row(i),mapper.reference_row(i))
        self.assertEqual(mapper.row(M.K.MAX_ID),mapper.reference_row(M.K.MAX_ID))
        for bad in (-1,2**32,True,1.0):
            with self.assertRaises(ValueError):mapper.row(bad)
        with self.assertRaisesRegex(ValueError,'noncommuting'):
            M.verify_mapper((pair[0],pair[0]),mapper,budget)
        with mock.patch.object(mapper,'payload',bytes(65535)):
            with self.assertRaisesRegex(ValueError,'geometry'):M.verify_mapper(pair,mapper,budget)

    def test_local_selection_first_success_and_exhaustion(self):
        records=[]
        with mock.patch.object(M.F,'matrix_rank',side_effect=[0,8]), \
             mock.patch.object(M,'trace') as trace,mock.patch.object(R,'inventory_inputs') as history:
            pair=M.G.choose_pair(NEUTRAL,('0',),(tuple(range(8)),),mock.Mock(),records)
        self.assertEqual([r['parameter'] for r in records],[1,2])
        self.assertEqual(pair[1],M.K.companion((7,)+NEUTRAL[1:]))
        trace.assert_not_called();history.assert_not_called()
        with mock.patch.object(M.F,'matrix_rank',return_value=0):
            records=[]
            self.assertIsNone(M.G.choose_pair(NEUTRAL,('0',),(tuple(range(8)),),mock.Mock(),records))
        self.assertEqual(len(records),254)
        self.assertNotIn(5,[r['parameter'] for r in records])

    def test_binary_factors_and_minor_roster(self):
        word='0'
        for _ in range(8):word=''.join('01' if c=='0' else '10' for c in word)
        self.assertEqual(set(M.WORDS),{word[i:i+4] for i in range(len(word)-3)})
        self.assertEqual(len(M.MINORS),495)
        self.assertEqual(len(set(M.MINORS)),495)
        self.assertEqual(M.MINORS,tuple(itertools.combinations(range(12),8)))

    def test_integer_trace_reference(self):
        def reference(width,root,kind):
            state=(int(root,16)^8*0x9e3779b97f4a7c15^width*0xbf58476d1ce4e5b9)%(2**64)
            if kind!='iid':state^=0x10fade
            threshold=(1/9 if kind=='burst' else .1 if kind=='iid' else .5)*2**53
            ids=[];skip=0
            for n in range(68608):
                if skip:skip-=1;continue
                state=(state+0x9e3779b97f4a7c15)%(2**64);v=state
                for shift,mult in ((30,0xbf58476d1ce4e5b9),(27,0x94d049bb133111eb)):
                    v=((v^(v>>shift))*mult)%(2**64)
                if ((v^(v>>31))>>11)<threshold:
                    if kind=='burst':skip=7
                else:
                    ids.append(2**32-1-2*n if kind=='adversarial' else 8+n if kind=='repair-only' else n)
                    if len(ids)==12:return ids
            raise AssertionError('trace exhausted')
        for b,s,r in itertools.product(M.WIDTHS,M.SCHEDULES,('0x0000000000000000','0x123456789abcdef0')):
            self.assertEqual(M.trace(b,r,s),reference(b,r,s))
        with mock.patch.object(M,'CANDIDATE_LIMIT',1):
            with self.assertRaisesRegex(ValueError,'candidate cap'):M.trace(64,'0x0000000000000000','iid')

    def test_each_fresh_cell_bound_and_invalid_ranks(self):
        rows=[dict(B=b,schedule=s,root=str(i),ranks=[8]*5)
              for b,s in itertools.product(M.WIDTHS,M.SCHEDULES) for i in range(100)]
        for offset in range(0,len(rows),100):
            rows[offset]['ranks']=[7,8,8,8,8]
            self.assertTrue(M.summarize_fresh(rows,100)['fresh_pass'])
            rows[offset+1]['ranks']=[7]*5
            result=M.summarize_fresh(rows,100)
            self.assertFalse(result['fresh_pass'])
            self.assertEqual(result['cells'][offset//100]['first_success'],[98,1,0,0,0,1])
            rows[offset]['ranks']=rows[offset+1]['ranks']=[8]*5
        for bad in (rows[:-1],rows+rows[:1],rows[:-1]+rows[:1]):
            with self.assertRaises(ValueError):M.summarize_fresh(bad,100)
        for bad in ([7]*4,[7,8,7,8,8],[6,8,8,8,8],[9]*5,[True]*5):
            rows[0]['ranks']=bad
            with self.assertRaisesRegex(ValueError,'nested'):M.summarize_fresh(rows,100)

    def test_collisions_and_exclusion_roster_without_scoring(self):
        with mock.patch.object(M.F,'digest',return_value='0'*64):
            with self.assertRaisesRegex(ValueError,'collision'):M.fresh_roots([])
        fake_hash=lambda s:format(int(s.decode().rsplit('/',1)[1])+1,'016x')+'0'*48
        with mock.patch.object(M.F,'digest',side_effect=fake_hash):
            roots=M.fresh_roots([])
            with self.assertRaisesRegex(ValueError,'collision'):M.fresh_roots(roots[:1])
        hashes=[]
        def recording_hash(value):
            hashes.append(value.decode());return hashlib.sha256(value).hexdigest()
        with mock.patch.object(M.F,'digest',side_effect=recording_hash):
            excluded=M.exclusions(dict(roots=['inventory-sentinel']))
        self.assertIn('inventory-sentinel',excluded)
        self.assertTrue(set(M.K.H.MAIN_ROOTS)<=set(excluded))
        self.assertEqual(len(hashes),4*512)
        self.assertEqual({s.split(':fresh/')[0] for s in hashes},{
            'wirehair.wh2.k2-thue-morse-r0','wirehair.wh2.k3-thue-morse-r0',
            'wirehair.wh2.k5-thue-morse-r0','wirehair.wh2.thue-morse-recovery-r0'})

    def pipeline(self,failure=None):
        pair=(M.K.companion(NEUTRAL),M.K.companion((9,)+NEUTRAL[1:]))
        class FakeMapper:
            def __init__(self,*args):self.cache={}
            def row(self,i):
                value=M.K.identity(8)[i%8];self.cache[i]=value;return value
        trace_calls=0
        def trace(b,root,schedule,mapper):
            nonlocal trace_calls
            trace_calls+=1
            fail=(failure=='hard' and trace_calls==1) or (failure=='fresh' and root in ('synthetic-0','synthetic-1','synthetic-2','synthetic-3','synthetic-4','synthetic-5'))
            return dict(B=b,root=root,schedule=schedule,ids=list(range(12)),ranks=[7,8,8,8,8] if fail else [8]*5)
        rank_calls=0
        def rank(rows):
            nonlocal rank_calls
            rank_calls+=1
            # Mapper verification is independently exercised with unrelated real matrices.
            # 4950 local minors, 14850 seam minors, then 42 historical prefixes.
            return 7 if (failure,rank_calls) in (('local',1),('seam',4951),('history',19801)) else 8
        with mock.patch.object(M,'fixed_feedback',return_value=NEUTRAL), \
             mock.patch.object(M.G,'choose_pair',return_value=pair) as choose, \
             mock.patch.object(M.K,'Mapper',FakeMapper),mock.patch.object(M,'verify_mapper',return_value={}), \
             mock.patch.object(M,'checked_rank',side_effect=rank), \
             mock.patch.object(R,'inventory_inputs',return_value=dict(prefixes=[dict(ids=list(range(8)),original_widths=[64])]*42,roots=[])), \
             mock.patch.object(M,'trace_result',side_effect=trace) as traces, \
             mock.patch.object(M,'fresh_roots',return_value=['synthetic-'+str(i) for i in range(512)]) as fresh:
            result=M.run_screen('synthetic')
        choose.assert_called_once()
        if failure=='local':
            self.assertEqual(result['outcome'],'INVALID');fresh.assert_not_called();traces.assert_not_called()
            return
        entered=failure not in ('seam','history','hard')
        self.assertEqual(result['outcome'],'PASS' if failure is None else 'FAIL',result.get('error'))
        self.assertEqual(traces.call_count,6216 if entered else 72)
        self.assertEqual(fresh.call_count,int(entered))
        self.assertEqual(len(result['fresh']),6144 if entered else 0)
        self.assertEqual(result['counts']['local_minors'],4950)
        self.assertEqual(result['counts']['seam_minors'],14850)
        self.assertEqual(result['counts']['history_prefixes'],42)
        self.assertEqual([r['ids'][0] for r in result['seams']],[(1<<e)-4 for e in range(3,32)]+[2**32-12])
        self.assertTrue(all(len(r['ids'])==12 for r in result['seams']))

    def test_full_synthetic_pipeline(self):self.pipeline()

    def test_structural_failure_stops_before_fresh_without_reselection(self):
        for failure in ('seam','history','hard','local'):
            with self.subTest(failure=failure):self.pipeline(failure)

    def test_fresh_failure_is_terminal_without_reselection(self):self.pipeline('fresh')

    def test_exhaustion_does_not_load_history_or_mapper(self):
        with mock.patch.object(M,'fixed_feedback',return_value=NEUTRAL), \
             mock.patch.object(M.G,'choose_pair',return_value=None), \
             mock.patch.object(R,'inventory_inputs') as history,mock.patch.object(M.K,'Mapper') as mapper:
            self.assertEqual(M.run_screen('synthetic')['outcome'],'EXHAUSTED')
        history.assert_not_called();mapper.assert_not_called()

    def test_claim_precedes_candidate(self):
        with mock.patch.object(R,'claimed_inputs',side_effect=ValueError('bad claim')), \
             mock.patch.object(M.resource,'setrlimit'),mock.patch.object(M,'run_screen') as run:
            with self.assertRaisesRegex(ValueError,'bad claim'):M.main(['--worker'])
            run.assert_not_called()
        for args in ([],['--run'],['--worker','extra']):
            with mock.patch.object(M,'run_screen') as run:
                with self.assertRaises(ValueError):M.main(args)
                run.assert_not_called()

    def test_actual_retained_projection_and_corruption(self):
        actual=R.inventory_inputs()
        self.assertEqual((len(actual['origins']),len(actual['prefixes']),len(actual['roots'])),(44,42,64))
        self.assertTrue(all(8<=len(p['ids'])<=12 and p['original_widths'] for p in actual['prefixes']))
        with mock.patch.object(R,'MANIFEST','0'*64):
            with self.assertRaisesRegex(ValueError,'manifest SHA'):R.inventory_inputs()
        raw=R.C.read_regular(R.INVENTORY/'raw.jsonl',64*1024**2)
        with self.assertRaises(ValueError):R.project(raw.replace(b'"ordinal":0,',b'"ordinal":1,',1))
        with mock.patch.dict(R.PROJECTIONS,origins=(44,'0'*64)):
            with self.assertRaisesRegex(ValueError,'projection origins'):R.project(raw)

    def test_receipt_inert_positive_claim_and_import_closure(self):
        receipt=dict(synthetic=True);history=dict(provenance=dict(synthetic=True))
        with mock.patch.object(R,'base_receipt',return_value=receipt), \
             mock.patch.object(R,'inventory_inputs',return_value=history),mock.patch.object(M,'run_screen') as run:
            actual=R.current_receipt()
            self.assertEqual(actual['inventory']['projection_sha256'],R.C.sha(R.C.canonical(history)))
            run.assert_not_called()
        claim=dict(protocol=R.C.PROTOCOL,receipt=receipt,receipt_sha256=R.C.sha(R.C.canonical(receipt)))
        with tempfile.TemporaryDirectory() as directory:
            R.C.write_new(Path(directory)/'CLAIM.json',R.C.canonical(claim))
            with mock.patch.object(R.C,'OUTPUT',Path(directory)),mock.patch.object(R,'current_receipt',return_value=receipt):
                self.assertEqual(R.claimed_inputs(),R.C.sha(R.C.canonical(claim)))
            with mock.patch.object(R.C,'OUTPUT',Path(directory)),mock.patch.object(R,'current_receipt',return_value={}):
                with self.assertRaises(ValueError):R.claimed_inputs()
        seen=set();paths=set()
        def visit(module):
            if id(module) in seen:return
            seen.add(id(module));name=getattr(module,'__file__',None)
            if name and R.C.ROOT in Path(name).resolve().parents:
                paths.add(str(Path(name).resolve().relative_to(R.C.ROOT)))
                for value in vars(module).values():
                    if isinstance(value,types.ModuleType):visit(value)
        visit(M)
        self.assertTrue(paths<=set(R.C.SOURCES),paths-set(R.C.SOURCES))


def load_tests(loader,tests,pattern):
    return unittest.TestSuite((tests,loader.loadTestsFromTestCase(OLD.FileTests),
        loader.loadTestsFromTestCase(OLD.CaptureTests),loader.loadTestsFromTestCase(OLD.PublicationTests)))


if __name__=='__main__':unittest.main()
