from pathlib import Path
import unittest
import Generate as G

ROOT = Path(__file__).resolve().parents[2]


class Test(unittest.TestCase):
    def test_exact_two_decoder_arithmetic_sites_only(self):
        raw = (ROOT/'codec/WirehairSmallCore.h').read_bytes()
        source = G.candidate(raw)
        original = raw.decode()
        for old,new in reversed(G.EDITS):
            self.assertEqual(source.count(new),1)
            source = source.replace(new,old)
        self.assertEqual(source,original)
        candidate = G.candidate(raw)
        self.assertEqual(candidate[:candidate.index('    Result Feed(')],
                         original[:original.index('    Result Feed(')])
        self.assertEqual(candidate[candidate.index('    Result Recover('):],
                         original[original.index('    Result Recover('):])
        self.assertLess(candidate.index('inverse = gf256_inv(row[pivot])'),
                        candidate.index('row[pivot] = 1;'))

    def test_fail_closed(self):
        raw = (ROOT/'codec/WirehairSmallCore.h').read_bytes()
        with self.assertRaises(ValueError): G.candidate(raw+b'\n')
        for source in ('','xx'):
            with self.assertRaises(ValueError): G.replace_once(source,'x','y')


if __name__ == '__main__':
    unittest.main()
