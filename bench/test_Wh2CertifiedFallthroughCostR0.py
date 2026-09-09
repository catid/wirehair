"""Neutral-only tests; never launch the timing worker."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

SPEC=importlib.util.spec_from_file_location('fallthrough_tested',
    Path(__file__).with_name('Wh2CertifiedFallthroughCostR0.py'))
C=importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(C)


def synthetic(ratio=1):
    records=[]
    for c in C.R.roster(C.CASES):
        duration=round(1000000*(ratio if c[5]==2 and C.R.SIDES[c[6]]^c[2] else 1))
        records.append(dict(coordinate=c,observation=dict(clocks=[0,0,0,duration,0,0])))
    return C.R.statistics(records,C.CASES)


class Tests(unittest.TestCase):
    def test_old_defaults(self):
        cfg=C.R.configuration(None)
        self.assertEqual(cfg.cases,C.R.CASES)
        self.assertIsNone(cfg.header_checker); self.assertIsNone(cfg.result_combiner)
        self.assertEqual(C.R.callback_count(),51840)
        self.assertEqual(len(list(C.R.roster())),51840)
        self.assertEqual(C.D.SETTINGS.cases,C.R.CASES)

    def test_extended_roster_and_phases(self):
        self.assertEqual(C.CASES[:20],C.R.CASES)
        self.assertEqual(len(C.CASES),len(set(C.CASES)))
        self.assertEqual(C.R.callback_count(C.CASES),82944)
        seen,phases={},{}
        for index,c in enumerate(C.R.roster(C.CASES)):
            i,r,o,w,m,pair,p,a,q=c
            self.assertEqual(i,index)
            self.assertEqual(a,C.R.PAIRS[pair][C.R.SIDES[p]^o])
            key=w,m,pair,o
            if p==0: seen.setdefault(key,[]).append(r)
            if p>=2 and p%2==0: phases.setdefault(key,[]).append((q*96//1000000)//2)
        self.assertEqual(len(seen),384)
        self.assertEqual(sum(k[2]<2 for k in seen),256)
        for key in seen:
            self.assertEqual(sorted(seen[key]),list(range(12)))
            self.assertEqual(sorted(phases[key]),sorted(list(range(48))*2))

    def test_small_rows_against_separate_oracles(self):
        k3=C.R.O.selected_rows()
        self.assertEqual(C.rows(3),k3[:11]+k3[12:])
        k5=C.R.sibling('fallthrough_k5_independent','Wh2K5PublicRecoveryR0.py')
        self.assertEqual(C.rows(5),tuple(k5.coefficient(i) for i in C.packet_ids(5)))
        for k in (3,5):
            self.assertEqual(C.rows(k)[:k],tuple(tuple(int(i==j) for j in range(k)) for i in range(k)))
            for b in (2,64,1280):
                arm=C.small_fixture(k,b)
                self.assertEqual(len(bytes.fromhex(arm['packets'])),(k+14)*b)
                self.assertGreaterEqual(arm['steps'],k)
                self.assertEqual(bytes.fromhex(arm['packets'])[:k*b],bytes((37*i+i//11)%256 for i in range(k*b)))

    def test_small_fixture_tampering(self):
        fixtures=[None]*20
        for c in C.CASES[20:]:
            _,k,b,_=c
            fixtures.append(dict(case=list(c),batch=128,
                source=bytes((37*i+i//11)%256 for i in range(k*b)).hex(),arms=[C.small_fixture(k,b)]*2))
        header=dict(fixtures=fixtures)
        with patch.object(C.R,'verify_header') as old:
            C.verify_header(header,0,'claim',[])
            self.assertEqual(len(old.call_args[0][0]['fixtures']),20)
            for slot in range(20,32):
                for key,value in (('case',[4,5,1280,2]),('batch',4),('source',''),('arms',[])):
                    bad=copy.deepcopy(header); bad['fixtures'][slot][key]=value
                    if bad==header: continue
                    with self.assertRaises(ValueError): C.verify_header(bad,0,'claim',[])
                for field in ('profile','packets','steps'):
                    bad=copy.deepcopy(header); bad['fixtures'][slot]['arms'][1][field]=0
                    with self.assertRaises(ValueError): C.verify_header(bad,0,'claim',[])
            with self.assertRaises(ValueError): C.verify_header(dict(fixtures=fixtures[:-1]),0,'claim',[])

    def test_retention_without_benefit_not_retained(self):
        same=synthetic()
        self.assertEqual(len(same['statistics']),384)
        result=C.combine([same,same])
        self.assertEqual(result['outcome'],'PASS')
        self.assertFalse(result['candidate_retained'])
        self.assertEqual(result['certified_four_cell_improvements'],[])

    def test_extended_raw_checks_all_cases_and_footer(self):
        zero=[0]*4
        previous=dict(clocks=[1,0,2,3,0,4],before=zero,after=zero)
        header=dict(prelude=previous,fixtures=[dict(arms=[dict(steps=c[1])]*2) for c in C.CASES])
        raw=[header]; addresses={n:[123]*n+[0]*(128-n) for n in (4,128)}
        for c in C.R.roster(C.CASES):
            ready=previous['clocks'][5]+1; cpu=previous['clocks'][4]
            target=ready+c[-1]; start=target+3; case=C.CASES[c[3]]; cycles=C.R.batch(case)
            k,metric=case[1],c[4]
            observed=dict(clocks=[start,cpu+2,start+1,start+1000001,cpu+1000002,start+1000002],before=zero,after=zero)
            raw.append(dict(type='record',coordinate=c,ready=ready,target=target,
                wait=[ready+1,cpu,target+2,cpu+1],observation=observed,
                counts=[0 if metric else cycles,0 if metric else cycles*(k+14),cycles if metric else 0,
                        cycles*k if metric else 0,cycles if metric else 0,cycles],
                addresses=addresses[cycles],address_count=cycles,complete=True,checked=True))
            previous=observed
        count=C.R.callback_count(C.CASES)
        raw.append(dict(type='footer',complete=True,records=count,work_ns=count*1000000))
        with patch.object(C,'verify_header') as checker:
            result=C.R.verify_rows(raw,'claim',1,[],C.PROTOCOL,C.CASES,checker)
            self.assertEqual(result['outcome'],'PASS')
            checker.assert_called_once_with(header,1,'claim',[],C.PROTOCOL)
            index=next(i for i,r in enumerate(raw[1:-1],1) if r['coordinate'][3]==31)
            saved=raw[index]
            raw[index]=dict(saved,counts=[0]*6)
            with self.assertRaisesRegex(ValueError,'every attempted API call'):
                C.R.verify_rows(raw,'claim',1,[],C.PROTOCOL,C.CASES,checker)
            raw[index]=saved; raw[-1]['records']=51840
            with self.assertRaisesRegex(ValueError,'terminal complete footer'):
                C.R.verify_rows(raw,'claim',1,[],C.PROTOCOL,C.CASES,checker)
            with self.assertRaisesRegex(ValueError,'whole raw cohort'):
                C.R.verify_rows(raw,'claim',1,[])

    def test_all_four_cells_and_controls_required(self):
        better=synthetic(.99); same=synthetic()
        self.assertTrue(C.combine([better,better])['candidate_retained'])
        self.assertFalse(C.combine([better,same])['candidate_retained'])
        for outcome in ('CONTROL_FAIL','REGRESSION','INCONCLUSIVE'):
            failed=copy.deepcopy(better); failed['outcome']=outcome
            result=C.combine([better,failed])
            self.assertEqual(result['outcome'],outcome)
            self.assertFalse(result['candidate_retained'])

    def test_work_body_unchanged(self):
        path='bench/Wh2AdmissionRegressionCostR0.cpp'
        old=C.R.command(['git','show','3af52cb:'+path]).decode()
        new=(C.ROOT/path).read_text()
        def work(text): return text[text.index('NOINLINE void RunWork('):text.index('void Prepare(')]
        self.assertEqual(hashlib.sha256(work(old).encode()).digest(),hashlib.sha256(work(new).encode()).digest())
        adapter=(C.ROOT/'bench/Wh2CertifiedFallthroughCostR0.cpp').read_text()
        self.assertIn('#include "Wh2CurrentPreservedDeferredCostR0.cpp"',adapter)
        self.assertNotIn('RunWork(',adapter)


if __name__=='__main__': unittest.main()
