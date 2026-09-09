"""Neutral tests only: no frozen K2 candidate or fresh observations."""
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
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


M=sibling('k2_tm_tested','Wh2K2ThueMorseR0.py'); R=M.R
OLD=sibling('k2_controller_tests','test_Wh2NoncommutingRadixRunR0.py'); OLD.M=R.C


class Tests(unittest.TestCase):
    def test_field_and_rank(self):
        M.F.init_field(); M.G.reference_rank([])
        self.assertEqual(M.G.REFERENCE[0],M.F.MUL); self.assertEqual(M.G.REFERENCE[1],M.F.INV)
        for i in range(256):
            raw=hashlib.sha256(('k2-neutral/'+str(i)).encode()).digest()
            rows=[raw[j:j+2] for j in range(0,12,2)]
            self.assertEqual(M.checked_rank(rows),M.G.reference_rank(rows))
        for r in ([],[(0,0)],[(1,0),(2,0)],[(1,0),(0,1)]):
            self.assertEqual(M.checked_rank(r),M.G.reference_rank(r))

    def test_pair_determinant_evidence(self):
        for i in range(256):
            raw=hashlib.sha256(('k2-pair-neutral/'+str(i)).encode()).digest()
            rows=[raw[:2],raw[2:4]]; result=M.checked_pair(rows)
            det=M.K.multiply_polynomial(rows[0][0],rows[1][1])^M.K.multiply_polynomial(rows[0][1],rows[1][0])
            self.assertEqual(result,dict(rank=M.G.reference_rank(rows),determinant=det))
        self.assertEqual(M.checked_pair([(1,0),(2,0)]),dict(rank=1,determinant=0))
        self.assertEqual(M.checked_pair([(0,0),(0,0)]),dict(rank=0,determinant=0))
        with mock.patch.object(M.K,'multiply_polynomial',return_value=0):
            with self.assertRaisesRegex(ValueError,'determinant'): M.checked_pair([(1,0),(0,1)])
        for rows in ([],[(1,0)],[(1,0),(0,1),(1,1)],[(1,),(0,)]):
            with self.assertRaises(ValueError): M.checked_pair(rows)

    def test_unrelated_dimension2_mapper(self):
        pair=(M.K.companion((5,7)),M.K.companion((9,11)))
        mapper=M.K.Mapper(pair,M.F.Budget()); product=M.K.identity(2)
        self.assertEqual(len(mapper.payload),7168)
        for i in range(2050):
            self.assertEqual(mapper.row(i),tuple(row[0] for row in product))
            product=M.F.matrix_multiply(product,pair[M.K.parity(i)])
        for bit in range(2,32):
            for offset in (-1,0,1): self.assertEqual(mapper.row((1<<bit)+offset),mapper.reference_row((1<<bit)+offset))
        self.assertEqual(mapper.row(M.K.MAX_ID),mapper.reference_row(M.K.MAX_ID))

    def test_local_first_success_and_exhaustion(self):
        records=[]
        with mock.patch.object(M.F,'matrix_rank',side_effect=[0,2]),mock.patch.object(M,'trace') as trace:
            pair=M.G.choose_pair((5,7),('0',),((0,1),),mock.Mock(),records)
        self.assertEqual([r['parameter'] for r in records],[1,2])
        self.assertEqual(pair[1],M.K.companion((7,7))); trace.assert_not_called()
        with mock.patch.object(M.F,'matrix_rank',return_value=0):
            records=[]; self.assertIsNone(M.G.choose_pair((5,7),('0',),((0,1),),mock.Mock(),records))
        self.assertEqual(len(records),254); self.assertNotIn(5,[r['parameter'] for r in records])

    def test_binary_factors_and_minor_count(self):
        text='0'
        for _ in range(8): text=''.join('01' if c=='0' else '10' for c in text)
        self.assertEqual(set(M.WORDS),{text[i:i+4] for i in range(len(text)-3)})
        self.assertEqual(len(M.MINORS),15)

    def test_trace_integer_reference(self):
        def reference(width,root,kind):
            state=(int(root,16)^2*0x9e3779b97f4a7c15^width*0xbf58476d1ce4e5b9)%(2**64)
            if kind!='iid': state^=0x10fade
            threshold=(1/9 if kind=='burst' else .1 if kind=='iid' else .5)*2**53
            ids=[]; skip=0
            for n in range(67072):
                if skip: skip-=1; continue
                state=(state+0x9e3779b97f4a7c15)%(2**64); v=state
                for shift,mult in ((30,0xbf58476d1ce4e5b9),(27,0x94d049bb133111eb)):
                    v=((v^(v>>shift))*mult)%(2**64)
                if ((v^(v>>31))>>11)<threshold:
                    if kind=='burst': skip=7
                else:
                    ids.append(2**32-1-2*n if kind=='adversarial' else 2+n if kind=='repair-only' else n)
                    if len(ids)==6: return ids
            raise AssertionError('trace exhausted')
        for b,s,r in itertools.product(M.WIDTHS,M.SCHEDULES,('0x0000000000000000','0x123456789abcdef0')):
            self.assertEqual(M.trace(b,r,s),reference(b,r,s))

    def test_fresh_cell_bound_and_collision(self):
        rows=[dict(B=b,schedule=s,root=str(i),ranks=[2]*5) for b,s in itertools.product(M.WIDTHS,M.SCHEDULES) for i in range(100)]
        for offset in range(0,len(rows),100):
            rows[offset]['ranks']=[1,2,2,2,2]; self.assertTrue(M.summarize_fresh(rows,100)['fresh_pass'])
            rows[offset+1]['ranks']=[1]*5; self.assertFalse(M.summarize_fresh(rows,100)['fresh_pass'])
            rows[offset]['ranks']=rows[offset+1]['ranks']=[2]*5
        for bad in (rows[:-1],rows+rows[:1],rows[:-1]+rows[:1]):
            with self.assertRaises(ValueError): M.summarize_fresh(bad,100)
        with mock.patch.object(M.F,'digest',return_value='0'*64):
            with self.assertRaisesRegex(ValueError,'collision'): M.fresh_roots([])
            with self.assertRaisesRegex(ValueError,'collision'): M.stride_pairs(257)
        with mock.patch.object(M.F,'digest',side_effect=lambda s:format(int(s.decode().rsplit('/',1)[1])+1,'016x')+'0'*48):
            for stride in M.STRIDES:
                pairs=M.stride_pairs(stride)
                self.assertEqual(pairs,[(i+1,i+1+stride) for i in range(512)])

    def pipeline(self,structural,failure=None):
        pair=(M.K.companion((1,0)),M.K.companion((1,1)))
        class FakeMapper:
            payload=bytes(7168)
            def __init__(self,*args):
                self.cache={}; self.low=[]; p=M.K.identity(2)
                for i in range(2049):
                    self.low.append(tuple(r[0] for r in p)); p=M.F.matrix_multiply(p,pair[M.K.parity(i)])
            def row(self,i):
                r=self.low[i] if i<2049 else M.K.identity(2)[i%2]; self.cache[i]=r; return r
        def trace(b,root,schedule,mapper):
            return dict(B=b,root=root,schedule=schedule,ids=list(range(6)),ranks=[2]*5 if structural else [1,2,2,2,2])
        rank_calls=0
        def rank(rows):
            nonlocal rank_calls
            rank_calls+=1
            # Two invertibility checks, 150 local minors, 450 seam minors,
            # then 13 historical prefixes and three legacy pairs.
            return 1 if (failure,rank_calls) in (('seam',153),('history',603),('legacy',616)) else 2
        with mock.patch.object(M,'fixed_feedback',return_value=(1,0)),mock.patch.object(M.G,'choose_pair',return_value=pair), \
             mock.patch.object(M.K,'Mapper',FakeMapper),mock.patch.object(M,'checked_rank',side_effect=rank), \
             mock.patch.object(M,'checked_pair',return_value=dict(rank=2,determinant=1)), \
             mock.patch.object(R,'inventory_inputs',return_value=dict(prefixes=[dict(ids=[0,1],original_widths=[64])]*13,roots=[])), \
             mock.patch.object(M,'trace_result',side_effect=trace) as traces, \
             mock.patch.object(M,'fresh_roots',return_value=['synthetic-'+str(i) for i in range(512)]) as fresh, \
             mock.patch.object(M,'stride_pairs',return_value=[(0,1)]*512) as strides:
            result=M.run_screen('synthetic')
        passed=structural and failure is None
        self.assertEqual(result['outcome'],'PASS' if passed else 'FAIL',result.get('error'))
        self.assertEqual(traces.call_count,6216 if passed else 72)
        self.assertEqual(fresh.call_count,int(passed)); self.assertEqual(strides.call_count,3 if passed else 0)
        self.assertEqual(len(result['fresh']),6144 if passed else 0)
        self.assertEqual(result['counts']['stride_pairs'],1536 if passed else 0)
        self.assertTrue(all('determinant' in pair for s in result['strides'] for pair in s['pairs']))

    def test_positive_full_synthetic_pipeline(self): self.pipeline(True)
    def test_hard_failure_stops_before_fresh(self): self.pipeline(False)

    def test_each_other_structural_failure_stops_before_fresh(self):
        for failure in ('seam','history','legacy'):
            with self.subTest(failure=failure): self.pipeline(True,failure)

    def test_exhaustion_does_not_load_history_or_mapper(self):
        with mock.patch.object(M,'fixed_feedback',return_value=(5,7)),mock.patch.object(M.G,'choose_pair',return_value=None), \
             mock.patch.object(R,'inventory_inputs') as history,mock.patch.object(M.K,'Mapper') as mapper:
            self.assertEqual(M.run_screen('synthetic')['outcome'],'EXHAUSTED')
        history.assert_not_called(); mapper.assert_not_called()

    def test_claim_precedes_candidate(self):
        with mock.patch.object(R,'claimed_inputs',side_effect=ValueError('bad claim')),mock.patch.object(M.resource,'setrlimit'), \
             mock.patch.object(M,'run_screen') as run:
            with self.assertRaisesRegex(ValueError,'bad claim'): M.main(['--worker'])
            run.assert_not_called()
        for args in ([],['--run'],['--worker','extra']):
            with mock.patch.object(M,'run_screen') as run:
                with self.assertRaises(ValueError): M.main(args)
                run.assert_not_called()

    def test_actual_retained_projection_and_corruption(self):
        actual=R.inventory_inputs()
        self.assertEqual((len(actual['origins']),len(actual['prefixes']),len(actual['roots'])),(56,13,64))
        with mock.patch.object(R,'MANIFEST','0'*64):
            with self.assertRaisesRegex(ValueError,'manifest SHA'): R.inventory_inputs()

    def test_positive_claim_and_import_closure(self):
        receipt=dict(synthetic=True)
        claim=dict(protocol=R.C.PROTOCOL,receipt=receipt,receipt_sha256=R.C.sha(R.C.canonical(receipt)))
        with tempfile.TemporaryDirectory() as directory:
            R.C.write_new(Path(directory)/'CLAIM.json',R.C.canonical(claim))
            with mock.patch.object(R.C,'OUTPUT',Path(directory)),mock.patch.object(R,'current_receipt',return_value=receipt):
                self.assertEqual(R.claimed_inputs(),R.C.sha(R.C.canonical(claim)))
            with mock.patch.object(R.C,'OUTPUT',Path(directory)),mock.patch.object(R,'current_receipt',return_value={}):
                with self.assertRaises(ValueError): R.claimed_inputs()
        seen=set(); paths=set()
        def visit(module):
            if id(module) in seen: return
            seen.add(id(module)); name=getattr(module,'__file__',None)
            if name and R.C.ROOT in Path(name).resolve().parents:
                paths.add(str(Path(name).resolve().relative_to(R.C.ROOT)))
                for v in vars(module).values():
                    if isinstance(v,types.ModuleType): visit(v)
        visit(M); self.assertTrue(paths<=set(R.C.SOURCES),paths-set(R.C.SOURCES))


def load_tests(loader,tests,pattern):
    return unittest.TestSuite((tests,loader.loadTestsFromTestCase(OLD.FileTests),
        loader.loadTestsFromTestCase(OLD.CaptureTests),loader.loadTestsFromTestCase(OLD.PublicationTests)))


if __name__=='__main__': unittest.main()
