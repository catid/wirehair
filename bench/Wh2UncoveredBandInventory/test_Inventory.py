"""Synthetic algebra, full reducer and API unwind checks; no native workload."""
import copy
import itertools
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import Inventory as I


class Test(unittest.TestCase):
    def test_complete_frozen_roster(self):
        rows = I.roster()
        self.assertEqual(I.KS,(7,9,12,16))
        self.assertEqual((len(rows),len(set(I.roots()))),(3096,64))
        self.assertEqual(sum(r['group']=='inventory' for r in rows),3072)
        self.assertEqual(sum(r['group']=='hard' for r in rows),24)
        self.assertEqual(len(I.neutral_roster()),16)
        self.assertFalse({r['b'] for r in rows}&{r['b'] for r in I.neutral_roster()})
        for row in rows:
            self.assertEqual(len(row['ids']),row['k']+4)
            self.assertEqual(len(set(row['ids'])),len(row['ids']))
            self.assertTrue(all(0<=value<=0xffffffff for value in row['ids']))

    def test_independent_trace(self):
        def reference(k,b,root,schedule):
            state = (int(root,16) ^ k*0x9e3779b97f4a7c15 ^ b*0xbf58476d1ce4e5b9)%(2**64)
            if schedule!='iid': state ^= 0x10fade
            threshold = (.5/(8-7*.5) if schedule=='burst' else .1 if schedule=='iid' else .5)*2**53
            out,skip,candidate = [],0,0
            while len(out)<k+4:
                if skip: skip -= 1
                else:
                    state=(state+0x9e3779b97f4a7c15)%(2**64); value=state
                    for shift,multiplier in ((30,0xbf58476d1ce4e5b9),(27,0x94d049bb133111eb)):
                        value=((value^(value>>shift))*multiplier)%(2**64)
                    word=(value^(value>>31))>>11
                    if word>=threshold:
                        out.append(0xffffffff-2*candidate if schedule=='adversarial' else
                                   k+candidate if schedule=='repair-only' else candidate)
                    elif schedule=='burst': skip=7
                candidate += 1
                self.assertLess(candidate,65536)
            return out
        for k,b,schedule,root in itertools.product(I.KS,I.WIDTHS,I.SCHEDULES,I.roots()):
            self.assertEqual(I.trace(k,b,root,schedule),reference(k,b,root,schedule))

    def test_exhaustive_field(self):
        for x in range(256):
            for y in range(256):
                a,b,value=x,y,0
                while b:
                    if b&1:value^=a
                    a=((a<<1)^(0x14d if a&128 else 0))&255; b>>=1
                self.assertEqual(I.products(x)[y],value)
            if x:self.assertEqual(I.products(x).count(1),1)

    def test_rank_and_packet_oracles(self):
        for k in I.KS:
            identity=[[int(i==j) for i in range(k)] for j in range(k)]
            self.assertEqual(I.rank(identity,k),k)
            self.assertEqual(I.rank(identity[:-1]+[identity[0]],k),k-1)
            self.assertEqual(I.first_rank([[0]*k]+identity,k),k+1)
            for b,tail in ((2,1),(17,17),(65,1)):
                data=I.source(k,b,tail)
                for index,row in enumerate(identity):
                    self.assertEqual(I.packet(row,data,b,index),data[index*b:(index+1)*b])
                row=list(range(1,k+1)); expected=[]
                for lane in range(b):
                    value=0
                    for index,coefficient in enumerate(row):
                        value ^= I.multiply(coefficient,data[index*b+lane] if index*b+lane<len(data) else 0)
                    expected.append(value)
                self.assertEqual(I.packet(row,data,b,0xffffffff),bytes(expected))
        for rows in ([[True,0]],[[256,0]],[[1]]):
            with self.assertRaises(ValueError):I.rank(rows,2)

    def test_summary_actual_endpoint_not_rank(self):
        records=[]
        for ordinal,row in enumerate(I.roster()):
            k=row['k']
            records.append(dict(ordinal=ordinal,case=row,arms=[
                dict(first=k+1 if k==9 else k,rank_first=k),dict(first=k,rank_first=k)]))
        result=I.summarize(records)
        self.assertEqual((result['recommended_k'],result['priority_order'][0]),(9,9))
        self.assertEqual(result['decoder_lag'],[774,0])
        self.assertEqual(result['totals'][1]['rank_failures'],[[0]*5,[0]*5])
        for key in ('candidate_tested','speed_claimed','holdout','all_K_claimed','promotion_claimed'):
            self.assertFalse(result[key])
        for row in records:row['arms'][0]['first']=row['case']['k']
        self.assertIsNone(I.summarize(records)['recommended_k'])

    def test_profile_identity(self):
        raw=I.struct.pack('<4sHHQQIB3s',b'WHV2',1,32,0x4b295bbb47f4f9c9,14,2,1,bytes(3))
        I.check_profile(raw.hex(),14,2)
        for bad in (raw[:8]+bytes(8)+raw[16:],raw[:-1]+b'\1',raw[:-1]):
            with self.assertRaises(ValueError):I.check_profile(bad.hex(),14,2)

    def test_encoder_error_unwinds(self):
        for exception in (False,True):
            api=I.Api.__new__(I.Api);api.arm='wh2';api.calls=[0]*6;freed=[]
            def create(source,message,b,profile,capacity,written,handle):
                handle._obj.value=42;written._obj.value=32
                raw=I.struct.pack('<4sHHQQIB3s',b'WHV2',1,32,0x4b295bbb47f4f9c9,message,b,0,bytes(3))
                I.T.memmove(profile,raw,32);return 0
            def encode(*args):
                if exception:raise RuntimeError('injected')
                return 2
            api.create=create;api.encode=encode;api.free=lambda handle:freed.append(handle.value)
            with self.assertRaises((ValueError,RuntimeError)):api.packets(bytes(14),2,[7])
            self.assertEqual(freed,[42]);self.assertEqual(api.calls,[1,1,0,0,0,1])

    def test_decoder_error_unwinds(self):
        api=I.Api.__new__(I.Api);api.arm='wh1';api.calls=[0]*6;freed=[]
        def create(unused,message,b,handle):handle._obj.value=43;return 0
        api.receiver=create;api.decode=lambda *args:2;api.free=lambda handle:freed.append(handle.value)
        with self.assertRaises(ValueError):api.receive('00'*32,bytes(14),2,[7],[bytes(2)])
        self.assertEqual(freed,[43]);self.assertEqual(api.calls,[0,0,1,1,0,1])

    def test_complete_synthetic_reducer(self):
        cases=I.neutral_roster()
        class Synthetic:
            def __init__(self,arm):self.arm=arm
            def packets(self,data,b,ids):
                k=(len(data)+b-1)//b
                profile=(I.struct.pack('<4sHHQQIB3s',b'WHV2',1,32,0x4b295bbb47f4f9c9,
                                       len(data),b,0,bytes(3)).hex() if self.arm=='wh2' else '00'*32)
                return profile,[data[(i%k)*b:(i%k+1)*b].ljust(b,b'\0')[
                    :len(data)-(k-1)*b if i==k-1 else b] for i in ids]
            def receive(self,profile,data,b,ids,packets):
                k=(len(data)+b-1)//b
                return dict(feed=[1]*(k-1)+[0],first=k,recovered=[I.A.sha(data)]*2)
        apis=[Synthetic(arm) for arm in ('wh2','wh1')]
        observed=I.observe_rows(apis,cases);records=[]
        I.evaluate(apis,cases,observed,records.append)
        raw=[dict(type='header',protocol=I.PROTOCOL,claim='test',library={},coefficients=observed)]+records+[
            dict(type='footer',complete=True,records=len(cases),calls=I.expected_calls(records,cases))]
        mutations=[lambda x:x[1].update(ordinal=1),
                   lambda x:x[1].update(ordinal=False),
                   lambda x:x[1]['case'].update(k=7.0),
                   lambda x:x[1]['arms'][0]['feed'].__setitem__(0,True),
                   lambda x:x[-1]['calls'][0].__setitem__(0,float(x[-1]['calls'][0][0])),
                   lambda x:x[1]['arms'][0].update(rank_first=0),
                   lambda x:x[-1]['calls'][0].__setitem__(5,0),
                   lambda x:x[1]['arms'][0]['packets'].__setitem__(0,'00'),
                   lambda x:x[0]['coefficients'][0]['ids'].__setitem__(0,42),
                   lambda x:x[1]['arms'][1].update(feed=[0])]
        with tempfile.TemporaryDirectory() as directory,mock.patch.object(I,'check_library'):
            path=Path(directory)/'raw'
            path.write_bytes(b''.join(I.A.canonical(row) for row in raw))
            result=I.verify(path,'test',neutral=True)
            self.assertEqual({k:v for k,v in result.items() if k!='parity_sha256'},
                             dict(cases=16,checked=True,scientific_launch=False,mode='native'))
            expected=I.hashlib.sha256()
            expected.update(I.A.canonical({k:raw[0][k] for k in ('type','protocol','coefficients')}))
            for row in raw[1:]: expected.update(I.A.canonical(row))
            self.assertEqual(result['parity_sha256'],expected.hexdigest())
            for first in (0,8):
                changed=copy.deepcopy(raw)
                arm=changed[1]['arms'][0]
                arm.update(first=first,feed=[1]*11 if first==0 else [1]*7+[0],
                           recovered=[] if first==0 else arm['recovered'])
                changed[-1]['calls']=I.expected_calls(changed[1:-1],cases)
                path.write_bytes(b''.join(I.A.canonical(row) for row in changed))
                with self.assertRaisesRegex(ValueError,'systematic neutral succeeds at K'):
                    I.verify(path,'test',neutral=True)
            for mutate in mutations:
                changed=copy.deepcopy(raw);mutate(changed)
                path.write_bytes(b''.join(I.A.canonical(row) for row in changed))
                with self.assertRaises(ValueError):I.verify(path,'test',neutral=True)
            path.write_bytes(b''.join(I.A.canonical(row) for row in raw)[:-1])
            with self.assertRaises(ValueError):I.verify(path,'test',neutral=True)


if __name__=='__main__':unittest.main()
