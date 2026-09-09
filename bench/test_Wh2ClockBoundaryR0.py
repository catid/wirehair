import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

SPEC = importlib.util.spec_from_file_location('clock_boundary', Path(__file__).with_name('Wh2ClockBoundaryR0.py'))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def fixture():
    cpuid = [[0,16,0x68747541,0x444d4163,0x69746e65], [1,0x00b00f81,100<<24,0,0],
             [0x80000000,0x80000028,0,0,0], [0x80000001,0,0,0,1<<27],
             [0x80000007,0,0,0,1<<8], [0x80000008,0,0,0,0], [0x80000021,4,0,0,0]]
    header = dict(type='header',protocol=M.PROTOCOL,claim='0'*64,neutral=True,samples=2,
                  iterations=M.ITERATIONS,seed=M.SEED,answer=0xd63028b03dfd593c,aux=50,
                  start=[500,0],cpuid=cpuid)
    rows = [[i,100+i*100,180+i*100,1000+i*100,1080+i*100,
             10000+i*1000,10020+i*1000,10240+i*1000,10260+i*1000,
             0xd63028b03dfd593c,0,0,0,0,0,0,0,0,50,50,50,50,11] for i in range(2)]
    return [header]+rows+[dict(type='footer',complete=True,records=2,codec_calls=0,end=[2000,500])]


class ReaderTest(unittest.TestCase):
    def analyze(self, values, failure=None):
        with tempfile.TemporaryDirectory(prefix='wh2-clock-reader-') as directory:
            path = Path(directory)/'raw.jsonl'
            path.write_text(''.join(json.dumps(v)+'\n' for v in values))
            return M.analyze(path,2,'0'*64,True,failure)

    def test_answer(self):
        self.assertEqual(M.answer(),0xd63028b03dfd593c)

    def test_valid(self):
        self.assertEqual(self.analyze(fixture())['outcome'],'NEUTRAL_PASS')

    def test_mutations(self):
        changes = [
            (0,'claim','1'*64), (0,'neutral',False), (0,'seed',0), (0,'iterations',1),
            (0,'answer',0), (0,'aux',51), (0,'samples',3),
            (1,0,1), (1,9,0), (1,22,10), (1,18,51), (1,19,51), (1,20,51), (1,21,51),
            (1,5,10030), (1,6,10250), (1,7,10270), (1,1,181), (1,3,1081),
            (1,10,1), (1,0,False), (2,1,0), (2,3,0), (2,5,0),
            (3,'complete',False), (3,'codec_calls',1), (3,'records',1), (3,'end',[0,0])]
        for index,key,value in changes:
            with self.subTest(index=index,key=key,value=value):
                data = fixture(); data[index][key] = value
                with self.assertRaises((ValueError,StopIteration)):
                    self.analyze(data)

    def test_identity(self):
        for leaf,word in [(0,2),(1,1),(1,2),(3,4),(4,4),(6,1)]:
            with self.subTest(leaf=leaf,word=word):
                data=fixture(); data[0]['cpuid'][leaf][word]=0
                with self.assertRaises(ValueError): self.analyze(data)

    def test_truncation_trailing(self):
        for values in [fixture()[:-1],fixture()[:-2],fixture()+[[]]]:
            with self.subTest(length=len(values)):
                with self.assertRaises((ValueError,StopIteration)): self.analyze(values)

    def test_counter_backwards_between(self):
        data=fixture(); data[1][11]=1
        with self.assertRaises(ValueError): self.analyze(data)

    def test_expected_result_failure(self):
        data=fixture(); data[2][9]^=1; data[3]['complete']=False
        self.assertEqual(self.analyze(data,'result')['outcome'],'EXPECTED_FAILURE')
        with self.assertRaises(ValueError): self.analyze(data)

    def test_expected_clock_failure(self):
        data=fixture(); data[2][22]=7
        for index in (2,4,8,21): data[2][index]=0
        data[3]['complete']=False
        self.assertEqual(self.analyze(data,'clock')['outcome'],'EXPECTED_FAILURE')
        with self.assertRaises(ValueError): self.analyze(data)
        broken=copy.deepcopy(data); broken[2][9]=0
        with self.assertRaises(ValueError): self.analyze(broken,'clock')

    def test_strict_json(self):
        for text in ['{"a":1,"a":2}','[NaN]','[Infinity]']:
            with self.assertRaises(ValueError): M.decode(text)

    def test_line_bound_and_unterminated(self):
        with tempfile.TemporaryDirectory(prefix='wh2-clock-lines-') as directory:
            path=Path(directory)/'bad'
            for value in ['[]',' '*8192+'\n']:
                path.write_text(value)
                with self.assertRaises(ValueError): list(M.records(path))

    def test_scientific_summary_and_missing_phase(self):
        values=fixture(); header=values[0]; header.update(neutral=False,samples=48)
        data=[header]
        for i in range(48):
            m=1000000+i*1000000+(2*i+1)*1000000//96
            compute=160000 if i==7 else 80000
            r=[i,m,m+compute+200,m,m+compute+20,
               (m-(60000 if i==8 else 10))*3,(m+10)*3,(m+10+compute)*3,(m+compute+30)*3,
               0xd63028b03dfd593c,0,0,0,0,0,0,0,0,50,50,50,50,11]
            data.append(r)
        data.append(dict(type='footer',complete=True,records=48,codec_calls=0,end=[60000000,60000000]))
        with tempfile.TemporaryDirectory(prefix='wh2-clock-summary-') as directory, mock.patch.object(M,'SAMPLES',48):
            path=Path(directory)/'raw.jsonl'
            path.write_text(''.join(json.dumps(v)+'\n' for v in data))
            result=M.analyze(path,48,'0'*64)
            self.assertEqual(result['phase_counts'],[1]*48)
            self.assertEqual([e['record'][0] for e in result['events']],[7,8])
            self.assertEqual([e['flags'] for e in result['events']],[[False,True,False],[True,False,False]])
            self.assertEqual(result['counts_ge_50us'],[1,48,0])
            self.assertEqual(result['counter_totals'],[0]*4)
            data[2][3]=data[2][3]//1000000*1000000+10416
            path.write_text(''.join(json.dumps(v)+'\n' for v in data))
            with self.assertRaisesRegex(ValueError,'48 passive'): M.analyze(path,48,'0'*64)


if __name__ == '__main__':
    unittest.main()
