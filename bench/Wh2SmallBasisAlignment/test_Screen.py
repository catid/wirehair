import io
import unittest
import Screen as S


def synthetic(control_failure=False, no_effect=False):
    yield ','.join(S.HEADER)+'\n'
    for r in S.roster():
        _,_,pair,_,_,side,carrier,offset = r
        duration = 100000
        if side and ((pair >= 4 and not no_effect) or (pair < 4 and control_failure)):
            duration = 125000
        pointers = [0x1000,0x2000,0x4000,0x1000000]
        row = r+[duration,pointers[0],pointers[1],0x100000+carrier*16384+offset,pointers[2],pointers[3]]
        yield ','.join(map(str,row))+'\n'


class Test(unittest.TestCase):
    def test_complete(self):
        r = S.analyze(synthetic())
        self.assertEqual(r['outcome'],'ALIGNMENT_EFFECT_DETECTED')
        self.assertEqual(r['rows'],34560)
        self.assertEqual(len(r['statistics']),160)
        self.assertEqual(sum(s['primary'] for s in r['statistics']),24)
        self.assertEqual(sum(s['control'] for s in r['statistics']),64)
        self.assertFalse(r['production_speed_claimed'])
        self.assertFalse(r['historical_slowdown_explained'])

    def test_controls_override_effect(self):
        r = S.analyze(synthetic(control_failure=True))
        self.assertEqual(r['outcome'],'CONTROL_FAIL')
        self.assertEqual(len(r['failed_controls']),64)

    def test_missing_effect(self):
        r = S.analyze(synthetic(no_effect=True))
        self.assertEqual(r['outcome'],'NO_UNIFORM_PRIMARY_EFFECT')
        self.assertEqual(len(r['missing_primary_effects']),24)

    def test_incomplete(self):
        with self.assertRaises(ValueError):
            S.analyze(io.StringIO(','.join(S.HEADER)+'\n'))

    def test_chronology_pointer_drift_and_overlap(self):
        for index,value in ((0,'1'),(8,'0'),(9,'0'),(9,'8192'),(10,'1088'),(11,'1048577'),(12,'1052672')):
            def changed():
                for i,line in enumerate(synthetic()):
                    if i == 1:
                        row = line.strip().split(','); row[index] = value
                        line = ','.join(row)+'\n'
                    yield line
            with self.assertRaises(ValueError):
                S.analyze(changed())

    def test_work_cap_includes_warmups(self):
        def oversized():
            for i,line in enumerate(synthetic()):
                if i:
                    row = line.strip().split(',')
                    if int(row[4]) < 2:
                        row[8] = '10000000000'
                    line = ','.join(row)+'\n'
                yield line
        with self.assertRaisesRegex(ValueError,'aggregate WORK cap'):
            S.analyze(oversized())

    def test_confidence(self):
        self.assertEqual(S.confidence([0]*12),dict(ratio=1.0,lower95=1.0,upper95=1.0))
        with self.assertRaises(ValueError):
            S.confidence([0]*11)


if __name__ == '__main__':
    unittest.main()
