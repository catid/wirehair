import io
import itertools
import math
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import Screen as S


def groups():
    return {(c,p,o):[math.log(0.95) if p >= 3 else 0.0]*12
            for c in range(210) for p in range(5) for o in range(2)}


def synthetic():
    yield ','.join(S.HEADER)+'\n'
    for row in S.roster():
        _,cell,pair,order,pos,arm = row
        ns = 95000 if pair >= 3 and S.SIDES[pos]^order else 100000
        # Retained warmups deliberately differ, so accidentally including them
        # would spoil controls. WH1's first-success count is independently kept.
        if pos < 2:
            ns *= pos+2
        steps = 0 if cell%3 == 0 else S.shape(cell)['k']+int(arm == 2)
        yield ','.join(map(str,row+[ns,steps,0x1000,0x1000000,0x10000000,0x10000100,0x20000000,0x30000000]))+'\n'


def changed(column, value, ordinal=1):
    for i,line in enumerate(synthetic()):
        if i == ordinal:
            row = line.strip().split(','); row[column] = value
            line = ','.join(row)+'\n'
        yield line


class Test(unittest.TestCase):
    def test_complete(self):
        r = S.analyze(synthetic())
        self.assertEqual(r['outcome'],'PASS')
        self.assertEqual((r['rows'],len(r['statistics']),len(r['primary'])),(453600,2100,4))
        self.assertEqual(len(r['fixed_addresses']),2520)
        self.assertEqual(len(r['first_success']),420)
        self.assertFalse(r['production_promotion_claimed'])
        for p in r['primary']:
            self.assertAlmostEqual(p['ratio'],0.95)
        for stat in r['statistics']:
            self.assertAlmostEqual(stat['ratio'],1.0 if stat['pair'] < 3 else 0.95)

    def test_roster(self):
        self.assertEqual(len(list(S.roster())),453600)
        self.assertEqual([S.shape(c)['route'] for c in (0,30,60,90,120,150,180)],list(range(7)))
        self.assertEqual([S.shape(c)['k'] for c in (0,30,60,90,120,150,180)],[3,3,5,8,3,16,3])
        self.assertEqual(sum(S.shape(c)['small'] for c in range(210)),150)
        for c in (-1,210,1.0):
            with self.assertRaises(ValueError): S.shape(c)

    def test_gate_precedence_and_wh1(self):
        g = groups()
        g[0,4,0] = [math.log(1.2)]*12
        r = S.summarize(g)
        self.assertEqual(r['outcome'],'FAIL')
        self.assertEqual(r['failed_required_wh1'],[[0,4,0]])
        self.assertEqual(r['candidate_not_proven_faster_than_WH1'],[[0,4,0]])
        g = groups()
        g[150,4,0] = [math.log(1.2)]*12  # certified WH1 deficit is reported only
        r = S.summarize(g)
        self.assertEqual(r['outcome'],'PASS')
        self.assertEqual(r['failed_required_wh1'],[])
        self.assertEqual(r['candidate_not_proven_faster_than_WH1'],[[150,4,0]])
        g[150,3,0] = [math.log(1.03)]*12  # certified control cannot be dropped
        self.assertEqual(S.summarize(g)['outcome'],'FAIL')
        g[209,2,1] = [math.log(1.03)]*12  # even WH1 AA overrides apparent wins
        self.assertEqual(S.summarize(g)['outcome'],'CONTROL_FAIL')

    def test_primary_missing_win(self):
        g = groups()
        for c in range(210):
            if S.shape(c)['small']: g[c,3,0] = [0.0]*12
        r = S.summarize(g)
        self.assertEqual(r['outcome'],'FAIL')
        self.assertEqual(sum(not p['passed'] for p in r['primary']),2)

    def test_encoder_not_primary_and_standalone_required(self):
        g = groups()
        for c in range(0,210,3):
            g[c,3,0] = g[c,3,1] = [math.log(1.01)]*12
        r = S.summarize(g)
        self.assertEqual(r['outcome'],'PASS')
        self.assertEqual({p['metric'] for p in r['primary']},{1,2})
        g[180,4,0] = [0.0]*12
        self.assertEqual(S.summarize(g)['outcome'],'FAIL')
        g = groups()
        # A route6 regression can erase the equal-fixture decoder mean win,
        # even with all per-cell retention bounds below1.02.
        for c in range(210):
            if S.shape(c)['small']:
                g[c,3,0] = [math.log(1.019 if c >= 180 else 0.999)]*12
        self.assertEqual(S.summarize(g)['outcome'],'FAIL')

    def test_covariance_preserved(self):
        g = groups()
        logs = [-0.03+0.02*(r%2) for r in range(12)]
        for c in range(210):
            if S.shape(c)['small']: g[c,3,0] = logs[:]
        r = S.summarize(g)
        expected = S.confidence(logs)
        for p in r['primary']:
            if p['order'] == 0:
                for key,value in expected.items(): self.assertAlmostEqual(p[key],value,places=14)

    def test_confidence_and_roster_reject(self):
        self.assertEqual(S.confidence([0]*12),dict(ratio=1.0,lower95=1.0,upper95=1.0))
        for logs in ([0]*11,[float('nan')]*12,[float('inf')]*12):
            with self.assertRaises(ValueError): S.confidence(logs)
        g = groups(); del g[209,4,1]
        with self.assertRaises(ValueError): S.summarize(g)

    def test_bad_rows_and_warmup_cap(self):
        for col,value in ((0,'1'),(6,'0'),(6,'-1'),(6,'nan'),(7,'1'),(8,'0'),(9,'16'),(11,str(0x10000000))):
            with self.assertRaises(ValueError): S.analyze(changed(col,value))
        with self.assertRaisesRegex(ValueError,'aggregate WORK cap'):
            S.analyze(changed(6,str(S.CAP-1)))
        with self.assertRaisesRegex(ValueError,'fixed lifecycle-cell'):
            S.analyze(changed(11,str(0x10000200),2))
        with self.assertRaisesRegex(ValueError,'complete numeric row'):
            S.analyze(itertools.islice(synthetic(),20))
        with self.assertRaises(ValueError): S.analyze(io.StringIO('wrong header\n'))

    def test_first_success_parity(self):
        # cell0 is encoder. cell1 begins after180 observations; alter its
        # second decoder sample while retaining the valid K..32 bound.
        with self.assertRaisesRegex(ValueError,'fixed arm first-success'):
            S.analyze(changed(7,'5',182))

    def test_strict_boundaries(self):
        g = groups(); g[0,0,0] = [math.log(1.02)]*12
        self.assertEqual(S.summarize(g)['outcome'],'CONTROL_FAIL')
        g = groups(); g[120,3,0] = [math.log(1.02)]*12
        self.assertEqual(S.summarize(g)['outcome'],'FAIL')

    def test_decisions(self):
        r = {k:{'outcome':'PASS'} for k in S.WORKERS}
        self.assertEqual(S.decisions(r)['outcome'],'PASS')
        r['normal']['outcome'] = 'FAIL'
        self.assertEqual(S.decisions(r)['outcome'],'FAIL')
        r['reverse']['outcome'] = 'CONTROL_FAIL'
        self.assertEqual(S.decisions(r)['outcome'],'CONTROL_FAIL')
        r['normal']['outcome'] = 'UNKNOWN'
        with self.assertRaises(ValueError): S.decisions(r)
        with self.assertRaises(ValueError): S.decisions({})

    def test_spent_namespace_before_subprocess(self):
        with tempfile.TemporaryDirectory() as temp, mock.patch.object(S,'OUTPUT',Path(temp)), \
                mock.patch.object(S,'inputs') as read_inputs:
            with self.assertRaisesRegex(ValueError,'namespace already spent'): S.run(Path(temp))
            read_inputs.assert_not_called()

    def test_sealed_manifests(self):
        with tempfile.TemporaryDirectory() as temp:
            member = Path(temp)/'binary'; member.write_bytes(b'qualified')
            seal = Path(temp)/'manifest'
            seal.write_text(S.digest(member)+'  '+str(member)+'\n')
            sha = S.digest(seal)
            self.assertEqual(S.checked_manifest(seal,sha),{str(seal):sha,str(member):S.digest(member)})
            with self.assertRaises(ValueError): S.checked_manifest(seal,'0'*64)
            member.write_bytes(b'replaced')
            with self.assertRaisesRegex(ValueError,'qualified artifact unchanged'): S.checked_manifest(seal,sha)
            seal.write_text(S.digest(member)+'  binary\n')
            self.assertEqual(S.checked_manifest(seal,S.digest(seal),Path(temp))[str(member)],S.digest(member))
            for content in ('','0'*64+'  relative\n',S.digest(member)+'  '+str(member)+'\n'+S.digest(member)+'  '+str(member)+'\n'):
                seal.write_text(content)
                with self.assertRaises(ValueError): S.checked_manifest(seal,S.digest(seal))


if __name__ == '__main__':
    unittest.main()
