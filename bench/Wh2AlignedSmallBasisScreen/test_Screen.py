import io
import unittest
import Screen as S


def synthetic(controlled, bad_control=False, missing_win=False, bad_retention=False):
    yield ','.join(S.CONTROLLED if controlled else S.NATURAL)+'\n'
    for row in S.roster(controlled):
        _,cell,pair,order,pos = row[:5]
        side = S.SIDES[pos]^order
        control = pair < (8 if controlled else 3)
        primary = (pair >= 8 and pair%4 != 0) if controlled else (pair == 3 and S.shape(cell)['k'] == 8 and
            S.shape(cell)['width'] == S.shape(cell)['tail'] == 1280 and cell%3 == 0)
        duration = 10000
        if side:
            if control:
                duration = 11000 if bad_control else 10000
            elif primary:
                duration = 10000 if missing_win else 8000
            else:
                duration = 11000 if bad_retention else 9900
        if controlled:
            _,_,_,_,_,_,arm,carrier,offset = row
            raw = 0x100000+carrier*16384+offset
            view = raw+(64-offset if arm and offset else 0)
            values = row+[duration,0x1000,0x2000,raw,view,0x4000,0x1000000]
        else:
            values = row+[duration,0x1000,0x1000000,0x10000000,0x20000000,0x30000000]
        yield ','.join(map(str,values))+'\n'


def changed(lines, column, value, ordinal=1):
    for i,line in enumerate(lines):
        if i == ordinal:
            row = line.strip().split(','); row[column] = value
            line = ','.join(row)+'\n'
        yield line


def trace_rows():
    yield ','.join(S.TRACE)+'\n'
    for k in (3,5,8):
        for width in (64,1280):
            for policy in range(2):
                for order in range(2):
                    for slot in range(2):
                        arm = slot^order
                        for allocation in range(3):
                            extra = 63 if arm and width == 1280 else 0
                            raw = 0x100000+slot*0x10000+allocation*0x4000+16
                            view = raw+((64-raw%64)%64 if allocation == 1 and extra else 0)
                            row = [k,width,policy,order,arm,allocation,(296,k*width+extra,48+8*k)[allocation],int(allocation == 1),
                                   raw,raw%64,raw%4096,view,view%64,view%4096]
                            yield ','.join(map(str,row))+'\n'


class Test(unittest.TestCase):
    def test_complete(self):
        for controlled,rows,intervals,aa,primary in ((True,41472,192,96,72),(False,298080,1380,828,4)):
            r = S.analyze(synthetic(controlled),controlled)
            self.assertEqual(r['outcome'],'PASS')
            self.assertEqual((r['rows'],len(r['statistics'])),(rows,intervals))
            self.assertEqual(sum(x['control'] for x in r['statistics']),aa)
            self.assertEqual(sum(x['primary'] for x in r['statistics']),primary)

    def test_controls_override(self):
        for controlled in (True,False):
            self.assertEqual(S.analyze(synthetic(controlled,bad_control=True),controlled)['outcome'],'CONTROL_FAIL')

    def test_missing_primary(self):
        for controlled in (True,False):
            self.assertEqual(S.analyze(synthetic(controlled,missing_win=True),controlled)['outcome'],'FAIL')

    def test_retention(self):
        for controlled in (True,False):
            self.assertEqual(S.analyze(synthetic(controlled,bad_retention=True),controlled)['outcome'],'FAIL')

    def test_bad_rows(self):
        for controlled,header,time,pointer in ((True,S.CONTROLLED,9,13),(False,S.NATURAL,6,8)):
            with self.assertRaises(ValueError):
                S.analyze(io.StringIO(','.join(header)+'\n'),controlled)
            for column,value in ((0,'1'),(time,'0'),(pointer,'0'),(pointer,'4096')):
                with self.assertRaises(ValueError):
                    S.analyze(changed(synthetic(controlled),column,value),controlled)
            with self.assertRaisesRegex(ValueError,'aggregate WORK cap'):
                S.analyze(changed(synthetic(controlled),time,'119999999999'),controlled)

    def test_view_and_pointer_drift(self):
        for column,value in ((12,'1048577'),(13,'1048640'),(14,'1088')):
            with self.assertRaises(ValueError):
                S.analyze(changed(synthetic(True),column,value),True)
        with self.assertRaisesRegex(ValueError,'fixed lifecycle-cell'):
            S.analyze(changed(synthetic(False),7,'8192',2),False)

    def test_confidence_and_decisions(self):
        self.assertEqual(S.confidence([0]*12),dict(ratio=1.0,lower95=1.0,upper95=1.0))
        for logs in ([0]*11,[float('nan')]*12):
            with self.assertRaises(ValueError):
                S.confidence(logs)
        with self.assertRaises(ValueError):
            S.decisions({})
        results = {k:{'outcome':'PASS'} for k in S.WORKERS}
        self.assertEqual(S.decisions(results)['outcome'],'PASS')
        results['controlled']['outcome'] = 'FAIL'
        self.assertEqual(S.decisions(results)['outcome'],'FAIL')
        results['natural']['outcome'] = 'CONTROL_FAIL'
        self.assertEqual(S.decisions(results)['outcome'],'CONTROL_FAIL')
        results['controlled']['outcome'] = 'UNKNOWN'
        with self.assertRaises(ValueError):
            S.decisions(results)

    def test_trace(self):
        r = S.trace_analyze(trace_rows())
        self.assertEqual((r['rows'],len(r['bases'])),(144,48))
        for column,value in ((6,'295'),(8,'0'),(11,'0')):
            with self.assertRaises(ValueError):
                S.trace_analyze(changed(trace_rows(),column,value))
        with self.assertRaises(ValueError):
            S.trace_analyze(io.StringIO(','.join(S.TRACE)+'\n'))


if __name__ == '__main__':
    unittest.main()
