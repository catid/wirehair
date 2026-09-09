#!/usr/bin/env python3
"""Neutral protocol/schema/decision tests. Never launch a scientific cohort."""
import copy
import importlib.util
from pathlib import Path
import unittest
from unittest import mock

SPEC = importlib.util.spec_from_file_location('admission_cost_tested',Path(__file__).with_name('Wh2AdmissionRegressionCostR0.py'))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.roster = list(M.roster())

    def records(self, ratio=1.0):
        result = []
        for c in self.roster:
            duration = round(1000000*(ratio if c[5]==2 and M.SIDES[c[6]]^c[2] else 1.0))
            result.append(dict(coordinate=c,observation=dict(clocks=[0,0,0,duration,0,0])))
        return result

    def test_frozen_roster(self):
        self.assertEqual(len(self.roster),51840)
        seen, phases = {}, {}
        for expected_index,c in enumerate(self.roster):
            index,r,o,w,m,pair,p,a,q = c
            self.assertEqual(index,expected_index)
            self.assertEqual(a,M.PAIRS[pair][M.SIDES[p]^o])
            key = (w,m,pair,o)
            if p==0:
                seen.setdefault(key,[]).append(r)
            if p>=2 and p%2==0:
                self.assertEqual(q,self.roster[index+1][-1])
                self.assertNotEqual(M.SIDES[p],M.SIDES[p+1])
                phase = (q*96//1000000)//2
                phases.setdefault(key,[]).append(phase)
        self.assertEqual(len(seen),240)
        for key,replicates in seen.items():
            self.assertEqual(sorted(replicates),list(range(12)))
            self.assertEqual(sorted(phases[key]),sorted(list(range(48))*2))

    def test_declared_cases(self):
        self.assertEqual(len(M.CASES),len(set(M.CASES)))
        self.assertEqual(len(M.CASES),20)
        self.assertEqual(sum(M.batch(c)==4 for c in M.CASES),2)
        self.assertEqual(M.CASES[:2],((0,2,2,2),(0,2,1280,2)))
        self.assertEqual(M.CASES[-2:],((3,6,2,2),(3,6,1280,2)))

    def test_equivalent_pass(self):
        result = M.statistics(self.records())
        self.assertEqual(result['outcome'],'PASS')
        self.assertEqual(len(result['statistics']),240)
        self.assertFalse(result['WH1_speed_qualified'])
        self.assertFalse(result['static_speed_qualified'])
        self.assertFalse(result['production_promotion_claimed'])

    def test_small_resolved_slowdown_fails(self):
        result = M.statistics(self.records(1.001))
        self.assertEqual(result['outcome'],'REGRESSION')
        self.assertEqual(len(result['resolved_regressions']),80)
        self.assertFalse(result['uncertain'])

    def test_improvement(self):
        result = M.statistics(self.records(.99))
        self.assertEqual(result['outcome'],'PASS')
        self.assertEqual(sum(r.get('resolved_improvement',False) for r in result['statistics']),80)

    def test_no_control_rescue(self):
        rows = self.records(.8)
        for row in rows:
            c = row['coordinate']
            if c[3:6]==[0,0,0] and M.SIDES[c[6]]^c[2]:
                row['observation']['clocks'][3]=1100000
        result = M.statistics(rows)
        self.assertEqual(result['outcome'],'CONTROL_FAIL')
        self.assertEqual(len(result['failed_controls']),2)
        self.assertFalse(result['resolved_regressions'])

    def test_uncertain_not_a_slowdown_claim(self):
        rows = self.records()
        for row in rows:
            c = row['coordinate']
            if c[5]==2 and M.SIDES[c[6]]^c[2]:
                row['observation']['clocks'][3]=1200000 if c[1]%2 else 850000
        result = M.statistics(rows)
        self.assertEqual(result['outcome'],'INCONCLUSIVE')
        self.assertEqual(len(result['uncertain']),80)
        self.assertFalse(result['resolved_regressions'])

    def test_incomplete_and_nonpositive_rejected(self):
        rows = self.records()
        with self.assertRaises(ValueError): M.statistics(rows[:-1])
        rows[2]['observation']['clocks'][3]=0
        with self.assertRaises(ValueError): M.statistics(rows)

    def test_corrupted_statistical_chronology_rejected(self):
        rows = self.records()
        rows[0] = copy.deepcopy(rows[0]); rows[0]['coordinate'][1] = True
        with self.assertRaises(ValueError): M.statistics(rows)

    def test_no_load_order_selection(self):
        for outcome in ('CONTROL_FAIL','REGRESSION','INCONCLUSIVE'):
            for order in (0,1):
                results = [dict(outcome='PASS'),dict(outcome='PASS')]
                results[order]['outcome'] = outcome
                combined = M.combine(results)
                self.assertEqual(combined['outcome'],outcome)
                self.assertFalse(combined['shared_preserved_path_screen_pass'])
        with self.assertRaises(ValueError): M.combine([dict(outcome='PASS')])

    def test_raw_ledger_chronology_and_failure_rejection(self):
        zero = [0]*4
        previous = dict(clocks=[1,0,2,3,0,4],before=zero,after=zero)
        header = dict(prelude=previous,fixtures=[dict(arms=[dict(steps=c[1])]*2) for c in M.CASES])
        rows = [header]
        addresses = {n:[123]*n+[0]*(128-n) for n in (4,128)}
        for c in self.roster:
            ready = previous['clocks'][5]+1; cpu = previous['clocks'][4]
            target = ready+c[-1]; start = target+3; cycles = M.batch(M.CASES[c[3]])
            k = M.CASES[c[3]][1]; metric = c[4]
            observed = dict(clocks=[start,cpu+2,start+1,start+1000001,cpu+1000002,start+1000002],before=zero,after=zero)
            rows.append(dict(type='record',coordinate=c,ready=ready,target=target,
                wait=[ready+1,cpu,target+2,cpu+1],observation=observed,
                counts=[0 if metric else cycles,0 if metric else cycles*(k+14),cycles if metric else 0,
                        cycles*k if metric else 0,cycles if metric else 0,cycles],
                addresses=addresses[cycles],address_count=cycles,complete=True,checked=True))
            previous = observed
        rows.append(dict(type='footer',complete=True,records=M.CALLBACKS,work_ns=M.CALLBACKS*1000000))
        with mock.patch.object(M,'verify_header'):
            self.assertEqual(M.verify_rows(rows,'0'*64,0,[])['outcome'],'PASS')
            saved = rows[1]
            for key,value in (('type','ignored'),('complete',False),('complete',1),('checked',False),
                              ('counts',[0]*6),('addresses',[0]*128),('address_count',1),
                              ('ready',0),('target',0),('wait',[0]*4)):
                rows[1] = copy.deepcopy(saved); rows[1][key] = value
                with self.assertRaises(ValueError): M.verify_rows(rows,'0'*64,0,[])
            rows[1] = copy.deepcopy(saved); rows[1]['observation']['clocks'][0]=0
            with self.assertRaises(ValueError): M.verify_rows(rows,'0'*64,0,[])
            rows[1] = saved; rows[-1]['complete']=False
            with self.assertRaises(ValueError): M.verify_rows(rows,'0'*64,0,[])

    def test_authenticated_metadata(self):
        meta = M.metadata()
        self.assertEqual(len(meta),2)
        for lib in meta:
            self.assertEqual(len(lib['exports']),53)
            self.assertEqual(len(lib['slots']),6)
            self.assertEqual(lib['context_bytes'],0x22810)
        header = M.bindings_header(meta)
        self.assertIn(b'const LibrarySpec library_specs[2]',header)
        self.assertEqual(header.count(b'libwirehair.so.2.0.0'),2)

    def test_binding_validation(self):
        meta = M.metadata()
        bindings = []
        for i,lib in enumerate(meta):
            base = (i+1)*0x1000000
            public = {s['name']:base+s['offset'] for s in lib['exports']}
            bindings.append(dict(base=base,context=base+lib['context'],
                                 exports=list(public.values()),slots=[public[s['name']] for s in lib['slots']],
                                 providers=[100,200,300,400,500],runtime_targets=list(range(1,38))))
        M.verify_bindings(bindings,meta)
        for field in ('context','exports','slots','providers','runtime_targets'):
            changed = copy.deepcopy(bindings)
            if field=='context': changed[1][field]+=1
            else: changed[1][field][0]+=1
            with self.assertRaises(ValueError): M.verify_bindings(changed,meta)

    def test_prior_fixture_coverage(self):
        prior = M.prior_records()
        for family,k,b,p in M.CASES:
            self.assertIn((M.FAMILIES[family],k,b,b,p),prior)


if __name__=='__main__':
    unittest.main()
