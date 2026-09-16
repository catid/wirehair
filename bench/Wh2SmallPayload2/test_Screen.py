import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import Screen as S


def synthetic(corrupt_control=False):
    yield ','.join(S.HEADER)+'\n'
    for r in S.roster():
        rep,cell,pair,order,pos,arm = r
        n = 200000 if arm == 2 else 80000 if arm == 1 and S.shape(cell)['width'] == 2 and cell%4 < 2 else 100000
        if corrupt_control and pair == 0 and (S.SIDES[pos]^order) == 1:
            n *= 2
        yield ','.join(map(str, r+[n]))+'\n'


class Test(unittest.TestCase):
    def test_complete_and_roster(self):
        result = S.analyze(synthetic())
        self.assertEqual(result['outcome'], 'PASS')
        self.assertEqual(len(result['statistics']), 960)
        self.assertEqual(result['candidate_not_proven_faster_than_WH1'], [])
        self.assertFalse(result['production_promotion_claimed'])

    def test_controls_override_wins(self):
        result = S.analyze(synthetic(True))
        self.assertEqual(result['outcome'], 'CONTROL_FAIL')
        self.assertEqual(len(result['failed_controls']), 192)

    def test_incomplete_or_shifted_or_zero(self):
        for text in (','.join(S.HEADER)+'\n', ','.join(S.HEADER)+'\n0,1,0,0,0,0,1\n',
                     ','.join(S.HEADER)+'\n0,0,0,0,0,0,0\n'):
            with self.assertRaises(ValueError):
                S.analyze(io.StringIO(text))

    def test_known_confidence(self):
        c = S.confidence([0]*12)
        self.assertEqual(c, dict(ratio=1.0, lower95=1.0, upper95=1.0))
        with self.assertRaises(ValueError):
            S.confidence([0]*11)

    def test_aggregate_work_bound_includes_warmups(self):
        def oversized():
            for line in synthetic():
                if line.startswith('rep,'):
                    yield line
                else:
                    row = line.strip().split(',')
                    if int(row[4]) < 2:
                        row[-1] = '10000000000'
                    yield ','.join(row)+'\n'
        with self.assertRaisesRegex(ValueError, 'aggregate WORK cap'):
            S.analyze(oversized())

    def test_failed_final_worker_cannot_replay_as_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'claim.json').write_text('{}')
            for mode in ('run','run-reverse'):
                (root/(mode+'.csv')).write_text('partial')
                (root/(mode+'.stderr')).write_text('')
                (root/(mode+'.json')).write_text(json.dumps(dict(exit=1, outcome='INVALID')))
            S.write(root/'complete.json', {p.name:S.digest(p) for p in root.iterdir()})
            with patch.object(S, 'analyze') as analysis:
                with self.assertRaisesRegex(ValueError, 'completed successfully'):
                    S.replay(root)
                analysis.assert_not_called()

    def test_spent_namespace_no_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'CMakeCache.txt').write_text('SANITIZE:BOOL=OFF\nPORTABLE:BOOL=OFF\nCOUNT:BOOL=OFF\n')
            (root/'libbaseline.so').write_text('baseline')
            (root/'libcandidate.so').write_text('candidate')
            with patch.object(S, 'OUTPUT', root), patch.object(S.subprocess, 'check_output', return_value=''), \
                    patch.object(S.subprocess, 'run') as run:
                with self.assertRaises(FileExistsError):
                    S.run(root)
                run.assert_not_called()

    def test_explicit_namespace_and_imported_sources(self):
        class StopBeforeLaunch(Exception):
            pass
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'CMakeCache.txt').write_text('SANITIZE:BOOL=OFF\nPORTABLE:BOOL=OFF\nCOUNT:BOOL=OFF\n')
            (root/'libbaseline.so').write_text('baseline')
            (root/'libcandidate.so').write_text('candidate')
            additional = root/'new-sources'
            additional.mkdir()
            (additional/'extra.py').write_text('# extra controller')
            output = root/'new-output'
            with patch.object(S.subprocess, 'check_output', return_value=''), \
                    patch.object(S, 'write', side_effect=StopBeforeLaunch) as write, \
                    patch.object(S.subprocess, 'run') as run:
                with self.assertRaises(StopBeforeLaunch):
                    S.run(root, output=output, here=additional, protocol='different-protocol')
                path, claim = write.call_args.args
                self.assertEqual(path, output/'claim.json')
                self.assertEqual(claim['protocol'], 'different-protocol')
                self.assertIn(str(additional/'extra.py'), claim['pins'])
                self.assertIn(str(S.HERE/'Screen.py'), claim['pins'])
                run.assert_not_called()


if __name__ == '__main__':
    unittest.main()
