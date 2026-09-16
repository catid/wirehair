"""Synthetic history checks and retired-wrapper preprocessing; no codec/timing."""
import io
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import Wh2K4RetainedHistory as H
import Wh2K4FreshCostBuildR1 as F


class HistoryTest(unittest.TestCase):
    def test_roster(self):
        rows = list(H.coordinates())
        self.assertEqual(len(rows),77760)
        self.assertEqual(rows[0],[0,0,0,0,0,0,0,0,10416])
        self.assertEqual(sum(r[6] < 2 for r in rows),8640)
        self.assertEqual(len({(r[1],r[2],r[3],r[4],r[5]) for r in rows}),4320)
        for start in range(0,len(rows),18):
            panel = rows[start:start+18]
            for pos in range(2,18,2):
                self.assertNotEqual(H.SIDES[panel[pos][6]],H.SIDES[panel[pos+1][6]])

    def test_confidence(self):
        c = H.confidence([0.0]*12)
        self.assertEqual((c['ratio'],c['lower95'],c['upper95']),(1.0,1.0,1.0))
        for values in ([0]*11,[math.nan]*12,[math.inf]*12):
            with self.assertRaises(ValueError): H.confidence(values)

    def test_compare_rejects_bad_schema_and_values(self):
        H.compare({'a':[0.5,1,True]}, {'a':[0.5+1e-15,1,True]})
        for other in ({'a':[0.5,True,True]}, {'a':[math.nan,1,True]}, {'a':[0.5,1]}, {'b':[]}):
            with self.assertRaises(ValueError): H.compare({'a':[0.5,1,True]},other)

    def test_json_rejects_duplicate_and_nonfinite(self):
        for raw in ('{"x":1,"x":2}','{"x":NaN}','{"x":Infinity}'):
            with self.assertRaises(ValueError): H.decode(raw)

    def test_seal_mismatch_stops_before_json(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); (root/'COMPLETE.json').write_bytes(b'not json')
            with self.assertRaisesRegex(ValueError,'fixed historical completion hash'):
                H.bundle('R22',root,'0'*64)

    def test_path_only_change_is_exact(self):
        names = 'bench/Wh2K4FreshCostBuildR1.py\nbench/Wh2K4SerializedCostR1.cpp\n'
        values = [names,b'claim R22\n',b'claim R23\n',b'claim R22\n',b'claim R23\n']
        with patch.object(H.subprocess,'check_output',side_effect=values):
            self.assertFalse(H.path_only_change('a802b66','6b2851e')['work_or_codec_change'])
        for values in ([b'other\n'],[names,b'claim R22\n',b'claim R23\nnew work\n']):
            with patch.object(H.subprocess,'check_output',side_effect=values), self.assertRaises(ValueError):
                H.path_only_change('a802b66','6b2851e')

    def test_bad_raw_coordinate_rejected(self):
        header = dict(type='header',protocol=H.PROTOCOL,claim='sealed',batch=128,fixtures=[{}]*6)
        record = dict(type='record',coordinate=[9]*9,complete=True,checked=True)
        stream = io.StringIO(H.canonical(header).decode()+H.canonical(record).decode())
        with patch.object(Path,'open',return_value=stream), self.assertRaisesRegex(ValueError,'raw chronology'):
            H.reconstruct(Path('/not-a-file'),'sealed',{})

    def test_spent_builder_and_qualifier_reject_before_io(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for mode in ('native','scalar','asan'):
                with patch.object(F,'_proof') as proof, patch.object(F.C,'build') as build:
                    with self.assertRaisesRegex(ValueError,'K4 R1 cost family retired'):
                        F.build(mode,root/'missing-proof',root/'missing-input',root/'new-output')
                    proof.assert_not_called(); build.assert_not_called()
                with patch.object(F.Q,'_neutral_root') as read:
                    with self.assertRaisesRegex(ValueError,'K4 R1 qualifier retired'):
                        F.Q.qualify(root/'missing-input',mode,root/'new-proof')
                    read.assert_not_called()
                self.assertEqual(list(root.iterdir()),[])

    def test_retired_wrapper_fails_before_including_codec(self):
        result = subprocess.run(['/usr/bin/c++','-E',str(H.ROOT/'bench/Wh2K4SerializedCostR1.cpp')],
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=10,check=False)
        self.assertNotEqual(result.returncode,0)
        self.assertIn(b'K4 R1 timing family retired',result.stderr)
        self.assertNotIn(b'Wh2SmallLifecycleWorkerR0.h',result.stdout)


if __name__ == '__main__': unittest.main()
