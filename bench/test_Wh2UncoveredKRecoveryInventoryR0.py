"""Synthetic inventory and cleanup checks; no scientific namespace or worker."""
import copy
import importlib.util
import itertools
from pathlib import Path
import tempfile
import types
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('inventory_tested',
    Path(__file__).with_name('Wh2UncoveredKRecoveryInventoryR0.py'))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class Tests(unittest.TestCase):
    def test_roster_exact_and_inert(self):
        r = M.roster()
        self.assertEqual(len(r),3096)
        self.assertEqual(len(set(M.roots())),64)
        self.assertEqual([sum(x['group']==g for x in r) for g in ('inventory','hard')],[3072,24])
        for k in M.KS:
            self.assertEqual(sum(x['group']=='inventory' and x['k']==k for x in r),768)
        for x in r:
            self.assertEqual(len(x['ids']),x['k']+4)
            self.assertEqual(len(set(x['ids'])),len(x['ids']))
            self.assertTrue(all(0<=i<=0xffffffff for i in x['ids']))

    def test_independent_trace_reference(self):
        # Different PRNG implementation and integer acceptance thresholds.
        def reference(k,b,root,schedule):
            state = (int(root,16) ^ (k*0x9e3779b97f4a7c15) ^ (b*0xbf58476d1ce4e5b9))%(2**64)
            if schedule!='iid': state ^= 0x10fade
            threshold = ((.5/(8-7*.5)) if schedule=='burst' else .1 if schedule=='iid' else .5)*2**53
            rows = []; skip = 0; candidate = 0
            while len(rows)<k+4:
                if skip: skip-=1
                else:
                    state=(state+0x9e3779b97f4a7c15)%(2**64)
                    v=state
                    for shift,multiplier in ((30,0xbf58476d1ce4e5b9),(27,0x94d049bb133111eb)):
                        v=((v^(v>>shift))*multiplier)%(2**64)
                    word=(v^(v>>31))>>11
                    if word>=threshold:
                        rows.append(0xffffffff-2*candidate if schedule=='adversarial' else k+candidate if schedule=='repair-only' else candidate)
                    elif schedule=='burst':skip=7
                candidate+=1
                self.assertLess(candidate,65536)
            return rows
        for k,b,schedule,root in itertools.product(M.KS,M.WIDTHS,M.SCHEDULES,M.roots()[:3]):
            self.assertEqual(M.trace(k,b,root,schedule),reference(k,b,root,schedule))

    def test_field_exhaustive_shift_reduce(self):
        for x in range(256):
            actual=M.products(x)
            for y in range(256):
                a,b,v=x,y,0
                while b:
                    if b&1:v^=a
                    a=((a<<1)^ (0x14d if a&128 else 0))&255;b>>=1
                self.assertEqual(actual[y],v)
            if x:self.assertEqual(actual.count(1),1)

    def test_rank_boundaries_and_decoder_lag(self):
        for k in M.KS:
            identity=[[int(i==j) for i in range(k)] for j in range(k)]
            self.assertEqual(M.rank(identity,k),k)
            self.assertEqual(M.rank(identity[:-1]+[identity[0]],k),k-1)
            self.assertEqual(M.first_rank([[0]*k]+identity,k),k+1)
        for bad in ([[256,0]],[[True,0]],[[1]]):
            with self.assertRaises(ValueError):M.rank(bad,2)

    def test_packet_full_tail_systematic_and_repair(self):
        for k,b,tail in itertools.product(M.KS,(2,64),(1,2)):
            data=M.source(k,b,tail)
            for i in range(k):
                row=[int(j==i) for j in range(k)]
                self.assertEqual(M.packet(row,data,b,i),data[i*b:(i+1)*b])
            row=list(range(1,k+1))
            expected=[]
            for lane in range(b):
                value=0
                for j,c in enumerate(row):
                    value ^= M.multiply(c,data[j*b+lane] if j*b+lane<len(data) else 0)
                expected.append(value)
            self.assertEqual(M.packet(row,data,b,0xffffffff),bytes(expected))

    def test_every_imported_project_helper_is_declared(self):
        seen=set(); paths=set()
        def visit(module):
            if id(module) in seen:return
            seen.add(id(module));name=getattr(module,'__file__',None)
            if name and M.ROOT in Path(name).resolve().parents:
                paths.add(str(Path(name).resolve().relative_to(M.ROOT)))
                for v in vars(module).values():
                    if isinstance(v,types.ModuleType):visit(v)
        visit(M)
        self.assertTrue(paths<=set(M.SOURCES))

    def test_bad_claim_is_rejected_before_loading(self):
        with patch.object(M.N,'Library') as load:
            with self.assertRaisesRegex(ValueError,'claim hex'):M.worker('bad')
            load.assert_not_called()

    def test_terminal_process_deadline_even_after_successful_exit(self):
        good=dict(returncode=0,error=None,elapsed_seconds=1.0)
        M.check_process(good)
        for elapsed in (0,150,151,float('nan'),float('inf'),True):
            with self.assertRaisesRegex(ValueError,'observer deadline'):
                M.check_process(dict(good,elapsed_seconds=elapsed))
        for code,error in ((1,None),(None,None),(False,None),(0,'observer error')):
            with self.assertRaisesRegex(ValueError,'worker completion'):
                M.check_process(dict(good,returncode=code,error=error))

    def test_custom_library_identity_rejects_before_dlopen(self):
        with patch.object(M.N.T,'CDLL') as load:
            with self.assertRaisesRegex(ValueError,'exact qualified DSO'):
                M.N.Library(0,((M.LIBRARY[0],'0'*64),))
            load.assert_not_called()

    def test_encoder_error_and_exception_free_owned_handle(self):
        for exception in (False,True):
            api=M.Api.__new__(M.Api);api.arm='wh2';api.calls=[0]*6;freed=[]
            def create(src,m,b,profile,capacity,written,handle):
                handle._obj.value=42;written._obj.value=32
                raw=M.struct.pack('<4sHHQQIB3s',b'WHV2',1,32,0x4b295bbb47f4f9c9,m,b,0,bytes(3))
                M.T.memmove(profile,raw,32);return 0
            def encode(*args):
                if exception:raise RuntimeError('injected')
                return 2
            api.create=create;api.encode=encode;api.free=lambda h:freed.append(h.value)
            with self.assertRaises((ValueError,RuntimeError)):api.packets(bytes(4),2,[2])
            self.assertEqual(freed,[42]);self.assertEqual(api.calls,[1,1,0,0,0,1])

    def test_decoder_error_frees_owned_handle(self):
        api=M.Api.__new__(M.Api);api.arm='wh1';api.calls=[0]*6;freed=[]
        def create(unused,m,b,h):h._obj.value=43;return 0
        api.receiver=create;api.decode=lambda *args:2;api.free=lambda h:freed.append(h.value)
        with self.assertRaisesRegex(ValueError,'decode status'):api.receive('00'*32,bytes(4),2,[2],[bytes(2)])
        self.assertEqual(freed,[43]);self.assertEqual(api.calls,[0,0,1,1,0,1])

    def test_summary_ranking_no_holdout_or_promotion(self):
        records=[]
        for ordinal,r in enumerate(M.roster()):
            k=r['k'];first=k+1 if k==5 else k
            records.append(dict(ordinal=ordinal,case=r,arms=[dict(first=first,rank_first=k),dict(first=k,rank_first=k)]))
        result=M.summarize(records)
        self.assertEqual(result['recommended_k'],5)
        self.assertEqual(result['priority_order'][0],5)
        self.assertFalse(result['holdout']);self.assertFalse(result['promotion_claimed'])
        self.assertEqual(result['decoder_lag'],[774,0])
        self.assertEqual(result['totals'][2]['rank_failures'],[[0]*5,[0]*5])
        for r in records:r['arms'][0]['first']=r['case']['k']
        self.assertIsNone(M.summarize(records)['recommended_k'])

    def test_profiles_reject_retired_and_wrong_shape(self):
        good=M.struct.pack('<4sHHQQIB3s',b'WHV2',1,32,0x4b295bbb47f4f9c9,8,2,1,bytes(3))
        M.check_profile(good.hex(),8,2)
        for bad in (good[:8]+bytes(8)+good[16:],good[:-1]+b'\1',good[:-1]):
            with self.assertRaises(ValueError):M.check_profile(bad.hex(),8,2)

    def test_full_reducer_path_and_tamper_failures(self):
        cases=[r for r in M.roster() if r['group']=='hard']
        class Synthetic:
            def __init__(self,arm):self.arm=arm
            def packets(self,data,b,ids):
                k=(len(data)+b-1)//b
                profile=(M.struct.pack('<4sHHQQIB3s',b'WHV2',1,32,0x4b295bbb47f4f9c9,len(data),b,0,bytes(3)).hex()
                         if self.arm=='wh2' else '00'*32)
                # Independent identity cycling, no native codec or field helper.
                return profile,[data[(i%k)*b:(i%k+1)*b].ljust(b,b'\0')[:len(data)-(k-1)*b if i==k-1 else b] for i in ids]
            def receive(self,profile,data,b,ids,packets):
                k=(len(data)+b-1)//b
                return dict(feed=[1]*(k-1)+[0],first=k,recovered=[M.A.sha(data)]*2)
        apis=[Synthetic(a) for a in ('wh2','wh1')]
        observed=M.observe_rows(apis,cases);records=[];M.evaluate(apis,cases,observed,records.append)
        elf=M.N.Elf(M.LIBRARY[0].read_bytes());base=0x1000000
        exports={name:dict(address=base+elf.symbol(name,2)[0],offset=elf.symbol(name,2)[0],size=elf.symbol(name,2)[1]) for name in elf.exports}
        context,size=elf.symbol('GF256Ctx',1)
        report=dict(path=str(M.LIBRARY[0]),sha256=M.LIBRARY[1],base=base,exports=exports,
            slots=[dict(name=n,offset=o,target=exports[n]['address']) for n,o in elf.slots],
            providers={n:base+1 for n in M.N.PROVIDERS},context=base+context,context_bytes=size,features=[1]*4)
        raw=[dict(type='header',protocol=M.PROTOCOL,claim='1'*64,library=report,coefficients=observed)]+records+[
            dict(type='footer',complete=True,records=len(cases),calls=M.expected_calls(records,cases))]
        mutations=[lambda x:x[1].update(ordinal=1),
            lambda x:x[1]['arms'][0].update(rank_first=0),
            lambda x:x[-1]['calls'][0].__setitem__(5,0),
            lambda x:x[0]['library']['slots'][0].update(target=0),
            lambda x:x[1]['arms'][0]['packets'].__setitem__(0,'00')]
        with tempfile.TemporaryDirectory() as folder,patch.object(M,'roster',return_value=cases),patch.object(M,'summarize',return_value={'verified':24}):
            path=Path(folder)/'synthetic.jsonl'
            path.write_bytes(b''.join(M.A.canonical(x) for x in raw))
            self.assertEqual(M.verify(path,'1'*64),{'verified':24})
            for change in mutations:
                altered=copy.deepcopy(raw);change(altered)
                path.write_bytes(b''.join(M.A.canonical(x) for x in altered))
                with self.assertRaises(ValueError):M.verify(path,'1'*64)
            path.write_bytes(b''.join(M.A.canonical(x) for x in raw)[:-1])
            with self.assertRaises(ValueError):M.verify(path,'1'*64)


if __name__=='__main__':unittest.main()
