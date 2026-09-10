"""Synthetic diagnostic tests. No production codec or scientific rerun."""
import copy
import hashlib
import math
import os
from pathlib import Path
import tempfile
import unittest

import Analyze as A


class Tests(unittest.TestCase):
    def test_exact_frozen_side_sequence(self):
        self.assertEqual(list(A.SIDES),[0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0])

    def fixture(self, order=0):
        header = dict(type='header',fixtures=[dict(case=[4,3,2,2])],
            bindings=[dict(base=4096,context=8192),dict(base=16384,context=32768)],
            prelude=dict(clocks=[1,0,2,3,1,4],before=[0]*4,after=[0]*4))
        rows = [header]; previous = header['prelude']; total = 0; logs = []
        for rep in range(12):
            contrasts = []
            for p in range(18):
                side = A.SIDES[p]^order; start = previous['clocks'][5]+20
                cpu = previous['clocks'][4]+10
                duration = 100+side*2+rep
                counters = [previous['after'][0]+2,0,0,previous['after'][3]]
                after = [counters[0]+1,0,0,counters[3]+(p==3)]
                row = dict(type='record',coordinate=[rep*18+p,rep,order,0,1,0,p,0,10],
                    ready=start-19,target=start-9,wait=[start-18,cpu-5,start-8,cpu-3],
                    observation=dict(clocks=[start,cpu,start+1,start+1+duration,cpu+duration+1,start+duration+3],
                                     before=counters,after=after),
                    complete=True,checked=True,address_count=2,addresses=[4096,4160,0,0])
                rows.append(row); previous=row['observation']; total+=duration
                if p>=2 and p%2:
                    left=A.SIDES[p-1]^order
                    values={left:100+left*2+rep,side:duration}
                    contrasts.append(math.log(values[1])-math.log(values[0]))
            logs.append(math.fsum(contrasts)/8)
        rows.append(dict(type='footer',complete=True,records=216,work_ns=total))
        stored=dict(outcome='CONTROL_FAIL',failed_controls=[[0,1,0,order]],
            statistics=[dict(case=[4,3,2,2],metric=1,comparison=0,order=order,replicate_logs=logs)])
        return rows,stored

    def test_all_records_and_preludes_preserved(self):
        rows,stored=self.fixture(); result=A.describe(iter(rows),stored,216)
        self.assertEqual(result['records'],216)
        self.assertEqual(len(result['cells']),1)
        cell=result['cells'][0]
        self.assertEqual(cell['original_statistics'],stored['statistics'][0])
        self.assertEqual(cell['same_address_vector_pairs'],108)
        for side in cell['sides']:
            self.assertEqual((side['records'],side['preludes']),(108,12))
            self.assertEqual(side['all_handle_mod64'],{0:216})
            self.assertEqual(side['first_addresses'],[4096])
        self.assertEqual(len(result['failed_control_panels']),1)
        self.assertEqual(sum(len(p['records']) for p in result['failed_control_panels'][0]['panels']),216)

    def test_clocks_and_counters_keep_separate_scopes(self):
        rows,_=self.fixture(); a=A.annotation(rows[1],rows[0]['prelude'])
        self.assertEqual(a['work_ns'],100)
        self.assertEqual(a['thread_bracket_ns'],101)
        self.assertEqual(a['capture_wall_ns'],103)
        self.assertEqual(a['wait_wall_ns'],10)
        self.assertEqual(a['wait_thread_ns'],2)
        self.assertEqual(a['wait_overshoot_ns'],1)
        self.assertEqual(a['capture_counters'],[1,0,0,0])
        self.assertEqual(a['gap_counters'],[2,0,0,0])

    def test_both_orders_keep_logical_ratio_and_all_failed_controls(self):
        results=[]
        for order in (0,1):
            rows,stored=self.fixture(order)
            result=A.describe(iter(rows),stored,216); results.append(result)
            self.assertTrue(all(x>0 for x in result['cells'][0]['original_statistics']['replicate_logs']))
            self.assertEqual(result['failed_control_panels'][0]['key'],[0,1,0,order])
            wrong=copy.deepcopy(stored); wrong['failed_controls']=[[9,1,0,order]]
            with self.assertRaisesRegex(ValueError,'every failed control'): A.describe(iter(rows),wrong,216)
        self.assertEqual(results[0]['cells'][0]['original_statistics']['replicate_logs'],
                         results[1]['cells'][0]['original_statistics']['replicate_logs'])

    def test_incomplete_changed_or_extra_records_reject(self):
        rows,stored=self.fixture()
        for broken in (rows[:-2], rows+[{}]):
            with self.assertRaises((ValueError,StopIteration)): A.describe(iter(broken),stored,216)
        for field,value in (('checked',False),('coordinate',[999]+rows[1]['coordinate'][1:])):
            broken=copy.deepcopy(rows); broken[1][field]=value
            with self.assertRaises(ValueError): A.describe(iter(broken),stored,216)
        changed=copy.deepcopy(stored); changed['statistics'][0]['replicate_logs'][0]+=1e-6
        with self.assertRaisesRegex(ValueError,'unchanged original paired statistics'):
            A.describe(iter(rows),changed,216)

    def test_address_vector_preserves_order_and_all_handles(self):
        rows,_=self.fixture(); first=A.annotation(rows[1],rows[0]['prelude'])
        changed=copy.deepcopy(rows[1]); changed['addresses']=[4160,4096,0,0]
        second=A.annotation(changed,rows[0]['prelude'])
        self.assertNotEqual(first['address_vector_sha256'],second['address_vector_sha256'])
        self.assertEqual(second['first_address'],4160)
        self.assertEqual(first['distinct_batch_addresses'],2)
        changed['addresses'][2]=1
        with self.assertRaises(ValueError): A.annotation(changed,rows[0]['prelude'])

    def test_record_counter_and_clock_reversals_reject(self):
        rows,_=self.fixture()
        for field in ('clocks','after'):
            changed=copy.deepcopy(rows[1]); changed['observation'][field][0]=0
            with self.assertRaises(ValueError): A.annotation(changed,rows[0]['prelude'])

    def test_exact_file_stream_and_json_reject_tampering(self):
        with tempfile.TemporaryDirectory(prefix='wh2-retained-diagnostic-test.') as folder:
            path=Path(folder)/'data'; raw=b'{"x":1}\n'
            fd=os.open(path,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
            try: os.write(fd,raw); os.fchmod(fd,0o400)
            finally: os.close(fd)
            pin=dict(path=str(path),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
            self.assertEqual(A.sealed_json(pin),{'x':1})
            with self.assertRaises(ValueError): list(A.sealed_lines(dict(pin,sha256='0'*64)))
            with self.assertRaises(ValueError): list(A.sealed_lines(pin,cap=1))
            link=Path(folder)/'link'; link.symlink_to(path)
            with self.assertRaises(OSError): list(A.sealed_lines(dict(pin,path=str(link))))
            os.chmod(path,0o600)
            with self.assertRaises(ValueError): A.sealed_json(pin)
            # Actual raw headers contain the complete 2.68-MB fixture roster.
            # A large synthetic line tests its transport, not codec arithmetic.
            large=Path(folder)/'large'; payload=b'"'+b'x'*(3*1024**2)+b'"\n'
            fd=os.open(large,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
            try:
                with os.fdopen(fd,'wb',closefd=False) as stream: stream.write(payload)
                os.fchmod(fd,0o400)
            finally: os.close(fd)
            large_pin=dict(path=str(large),bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest())
            self.assertEqual(b''.join(A.sealed_lines(large_pin)),payload)
        for raw in ('{"x":1,"x":2}','{"x":NaN}','{"x":Infinity}'):
            with self.assertRaises(ValueError): A.decode(raw)


if __name__=='__main__': unittest.main()
