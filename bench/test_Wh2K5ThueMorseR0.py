"""Neutral K5 arithmetic/controller tests; no frozen candidate selection."""
import contextlib
import hashlib
import importlib.util
import io
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


M=sibling('k5_tested','Wh2K5ThueMorseR0.py')
W=M.R
OLD=sibling('k5_old_controller_tests','test_Wh2NoncommutingRadixRunR0.py')
OLD.M=W.C


class Tests(unittest.TestCase):
    def test_independent_field_and_rank(self):
        M.F.init_field();M.reference_rank([])
        self.assertEqual(M.REFERENCE[0],M.F.MUL)
        self.assertEqual(M.REFERENCE[1],M.F.INV)
        for n in (2,4,5):
            identity=M.K.identity(n)
            for rank in range(n+1):
                rows=list(identity[:rank])+[(0,)*n]*3
                self.assertEqual(M.checked_rank(rows),rank)
        for i in range(128):
            raw=hashlib.sha256(('k5-neutral-matrix/'+str(i)).encode()).digest()
            self.assertEqual(M.reference_rank([raw[5*j:5*j+5] for j in range(5)]),
                             M.F.matrix_rank([raw[5*j:5*j+5] for j in range(5)]))
        for bad in ([[256,0]],[[True,0]],[[1],[0,1]]):
            with self.assertRaises(ValueError):M.reference_rank(bad)

    def test_selector_is_local_first_success_or_exhaustion(self):
        records=[]
        with mock.patch.object(M.F,'matrix_rank',side_effect=[0,3]), \
             mock.patch.object(M,'trace') as trace, mock.patch.object(W,'inventory_inputs') as history:
            pair=M.choose_pair((3,5,11),('0',),((0,1,2),),mock.Mock(),records)
        self.assertEqual([r['parameter'] for r in records],[1,2])
        self.assertIsNotNone(records[0]['first_failure']);self.assertIsNone(records[1]['first_failure'])
        self.assertEqual(pair[1],M.K.companion((1,5,11)))
        trace.assert_not_called();history.assert_not_called()
        with mock.patch.object(M.F,'matrix_rank',return_value=0):
            records=[]
            self.assertIsNone(M.choose_pair((3,5,11),('0',),((0,1,2),),mock.Mock(),records))
        self.assertEqual(len(records),254);self.assertNotIn(3,[r['parameter'] for r in records])

    def test_unrelated_dimension5_mapper_sequential_seams(self):
        pair=(M.K.companion((3,5,11,13,17)),M.K.companion((19,5,11,13,17)))
        mapper=M.K.Mapper(pair,M.F.Budget());product=M.K.identity(5)
        self.assertEqual(len(mapper.payload),29440)
        for i in range(2051):
            self.assertEqual(mapper.row(i),tuple(row[0] for row in product))
            product=M.F.matrix_multiply(product,pair[M.K.parity(i)])
        for bit in range(2,32):
            for offset in (-1,0,1):self.assertEqual(mapper.row((1<<bit)+offset),mapper.reference_row((1<<bit)+offset))
        self.assertEqual(mapper.row(M.K.MAX_ID),mapper.reference_row(M.K.MAX_ID))
        self.assertEqual(tuple(mapper.row(i) for i in range(5)),M.K.identity(5))

    def test_integer_rng_trace_reference(self):
        def reference(b,root,kind):
            state=(int(root,16)^5*0x9e3779b97f4a7c15^b*0xbf58476d1ce4e5b9)%(2**64)
            if kind!='iid':state^=0x10fade
            threshold=(1/9 if kind=='burst' else .1 if kind=='iid' else .5)*2**53
            ids=[];skip=0
            for n in range(67840):
                if skip:skip-=1;continue
                state=(state+0x9e3779b97f4a7c15)%(2**64);v=state
                for shift,mult in ((30,0xbf58476d1ce4e5b9),(27,0x94d049bb133111eb)):
                    v=((v^(v>>shift))*mult)%(2**64)
                if ((v^(v>>31))>>11)<threshold:
                    if kind=='burst':skip=7
                else:
                    ids.append(2**32-1-2*n if kind=='adversarial' else 5+n if kind=='repair-only' else n)
                    if len(ids)==9:return ids
            raise AssertionError('trace bound')
        with mock.patch.object(M,'choose_pair') as choose:
            for b,kind,root in itertools.product(M.WIDTHS,M.SCHEDULES,('0x0000000000000000','0x123456789abcdef0')):
                self.assertEqual(M.trace(b,root,kind),reference(b,root,kind))
        choose.assert_not_called()

    def test_fresh_gate_each_cell_exact_boundary_and_duplicates(self):
        rows=[dict(B=b,schedule=s,root=str(i),ranks=[5]*5) for b,s in itertools.product(M.WIDTHS,M.SCHEDULES) for i in range(100)]
        for offset in range(0,len(rows),100):
            rows[offset]['ranks']=[4,5,5,5,5]
            self.assertTrue(M.summarize_fresh(rows,100)['fresh_pass'])
            rows[offset+1]['ranks']=[4]*5
            result=M.summarize_fresh(rows,100);self.assertFalse(result['fresh_pass'])
            self.assertEqual(result['cells'][offset//100]['first_success'],[98,1,0,0,0,1])
            rows[offset]['ranks']=rows[offset+1]['ranks']=[5]*5
        for bad in (rows[:-1],rows+rows[:1],rows[:-1]+rows[:1]):
            with self.assertRaises(ValueError):M.summarize_fresh(bad,100)
        roots=M.fresh_roots([])
        with self.assertRaisesRegex(ValueError,'collision'):M.fresh_roots(roots[:1])

    def test_structural_failure_stops_before_fresh(self):
        pair=(M.K.companion((1,0,0,0,0)),)*2
        class FakeMapper:
            payload=bytes(29440)
            def __init__(self,*args):self.cache={}
            def row(self,i):
                value=M.K.identity(5)[i%5];self.cache[i]=value;return value
        def trace(b,root,schedule,mapper):return dict(B=b,root=root,schedule=schedule,ids=list(range(9)),ranks=[4,5,5,5,5])
        with mock.patch.object(M,'fixed_feedback',return_value=(1,0,0,0,0)), \
             mock.patch.object(M,'choose_pair',return_value=pair),mock.patch.object(M.K,'Mapper',FakeMapper), \
             mock.patch.object(M,'checked_rank',return_value=5),mock.patch.object(M,'trace_result',side_effect=trace) as traces, \
             mock.patch.object(W,'inventory_inputs',return_value={'prefixes':[list(range(5))]*54,'roots':[]}), \
             mock.patch.object(M,'fresh_roots') as fresh:
            result=M.run_screen('synthetic')
        self.assertEqual(result['outcome'],'FAIL',result.get('error'))
        self.assertFalse(result['summary']['fresh_entered']);self.assertEqual(result['fresh'],[])
        self.assertEqual(traces.call_count,72);fresh.assert_not_called()

    def test_positive_pipeline_enters_exact_fresh_roster_once(self):
        pair=(M.K.companion((1,0,0,0,0)),)*2
        class FakeMapper:
            payload=bytes(29440)
            def __init__(self,*args):self.cache={}
            def row(self,i):
                value=M.K.identity(5)[i%5];self.cache[i]=value;return value
        def trace(b,root,schedule,mapper):return dict(B=b,root=root,schedule=schedule,ids=list(range(9)),ranks=[5]*5)
        with mock.patch.object(M,'fixed_feedback',return_value=(1,0,0,0,0)), \
             mock.patch.object(M,'choose_pair',return_value=pair) as choose,mock.patch.object(M.K,'Mapper',FakeMapper), \
             mock.patch.object(M,'checked_rank',return_value=5),mock.patch.object(M,'trace_result',side_effect=trace) as traces, \
             mock.patch.object(W,'inventory_inputs',return_value={'prefixes':[list(range(5))]*54,'roots':[]}), \
             mock.patch.object(M,'fresh_roots',return_value=['synthetic-'+str(i) for i in range(512)]) as fresh:
            result=M.run_screen('synthetic')
        self.assertEqual(result['outcome'],'PASS',result.get('error'))
        self.assertTrue(result['summary']['fresh_entered']);self.assertEqual(len(result['fresh']),6144)
        self.assertEqual(result['summary']['failures'],[0]*5)
        self.assertEqual(traces.call_count,72+6144);fresh.assert_called_once();choose.assert_called_once()

    def test_binary_factor_roster_without_equations(self):
        word='0'
        for _ in range(8):word=''.join('01' if c=='0' else '10' for c in word)
        self.assertEqual({word[i:i+4] for i in range(len(word)-3)},set(M.WORDS))
        self.assertEqual(len(M.MINORS),126)

    def test_exhausted_never_loads_history_or_maps(self):
        with mock.patch.object(M,'fixed_feedback',return_value=(3,5,11)),mock.patch.object(M,'choose_pair',return_value=None), \
             mock.patch.object(W,'inventory_inputs') as history,mock.patch.object(M.K,'Mapper') as mapper:
            result=M.run_screen('synthetic')
        self.assertEqual(result['outcome'],'EXHAUSTED');history.assert_not_called();mapper.assert_not_called()

    def test_claim_authentication_precedes_selection(self):
        with mock.patch.object(W,'claimed_inputs',side_effect=ValueError('wrong claim')), \
             mock.patch.object(M.resource,'setrlimit'),mock.patch.object(M,'run_screen') as run:
            with self.assertRaisesRegex(ValueError,'wrong claim'):M.main(['--worker'])
            run.assert_not_called()
        for argv in ([],['--run'],['--worker','extra']):
            with mock.patch.object(M,'run_screen') as run:
                with self.assertRaisesRegex(ValueError,'usage'):M.main(argv)
                run.assert_not_called()

    def test_positive_claim_binding_and_mismatch(self):
        receipt={'sentinel':'synthetic'}
        claim=dict(protocol=W.C.PROTOCOL,receipt=receipt,receipt_sha256=W.C.sha(W.C.canonical(receipt)))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'CLAIM.json';W.C.write_new(path,W.C.canonical(claim))
            with mock.patch.object(W.C,'OUTPUT',Path(directory)),mock.patch.object(W,'current_receipt',return_value=receipt):
                self.assertEqual(W.claimed_inputs(),W.C.sha(W.C.canonical(claim)))
            with mock.patch.object(W.C,'OUTPUT',Path(directory)),mock.patch.object(W,'current_receipt',return_value={}):
                with self.assertRaisesRegex(ValueError,'authenticated'):W.claimed_inputs()

    def test_receipt_is_inert_and_all_imported_sources_declared(self):
        history={'provenance':{'synthetic':True}}
        with mock.patch.object(W,'_base_receipt',return_value={}),mock.patch.object(W,'inventory_inputs',return_value=history), \
             mock.patch.object(M,'run_screen') as run:
            self.assertEqual(W.current_receipt(),{'inventory':dict(provenance=history['provenance'],projection_sha256=W.C.sha(W.C.canonical(history)))})
            run.assert_not_called()
        seen=set();paths=set()
        def visit(module):
            if id(module) in seen:return
            seen.add(id(module));name=getattr(module,'__file__',None)
            if name and W.C.ROOT in Path(name).resolve().parents:
                paths.add(str(Path(name).resolve().relative_to(W.C.ROOT)))
                for v in vars(module).values():
                    if isinstance(v,types.ModuleType):visit(v)
        visit(M);self.assertTrue(paths<=set(W.C.SOURCES))

    def test_project_synthetic_origins_order_and_tamper(self):
        rows=[dict(type='case',ordinal=7,case=dict(k=5,b=64,root='synthetic',ids=list(range(9))),
            arms=[dict(first=7),dict(first=5)])]
        origins=[dict(ordinal=7,arm=0,b=64,overhead=oh,ids=list(range(5+oh))) for oh in (0,1)]
        expected=dict(origins=origins,prefixes=[list(range(5)),list(range(6))],roots=['synthetic'])
        pins={k:(len(v),W.C.sha(W.C.canonical(v))) for k,v in expected.items()}
        with mock.patch.object(W,'PROJECTIONS',pins):
            raw=b''.join(W.C.canonical(r)+b'\n' for r in rows)
            self.assertEqual(W.project(raw),expected)
            with self.assertRaises(ValueError):W.project(raw.replace(b'"first":7',b'"first":6'))


def load_tests(loader,tests,pattern):
    return unittest.TestSuite((tests,loader.loadTestsFromTestCase(OLD.FileTests),
        loader.loadTestsFromTestCase(OLD.CaptureTests),loader.loadTestsFromTestCase(OLD.PublicationTests)))


if __name__=='__main__':unittest.main()
