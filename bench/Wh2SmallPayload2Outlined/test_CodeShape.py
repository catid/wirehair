import unittest
import CodeShape as S


class Test(unittest.TestCase):
    def fixture(self):
        baseline, candidate = bytearray(64), bytearray(80)
        calls = []
        for i, owner in enumerate([S.ENCODER]*3 + ['wirehair_small_encode']):
            pc = i*8
            baseline[pc:pc+5] = b'\xe8' + (200-pc-5).to_bytes(4, 'little', signed=True)
            candidate[pc:pc+5] = b'\xe8' + (64-pc-5).to_bytes(4, 'little', signed=True)
            calls.append((pc,64,owner))
        return baseline, candidate, calls

    def compare(self, baseline, candidate, calls):
        return S.compare(0,baseline,candidate,calls,200,[(0,24)])

    def test_exact_redirects(self):
        self.assertEqual(self.compare(*self.fixture()), [])

    def test_other_difference_is_reported(self):
        b,c,calls = self.fixture()
        c[40] = 1
        self.assertEqual(self.compare(b,c,calls), [40])

    def test_protected_difference_is_rejected(self):
        b,c,calls = self.fixture()
        c[7] = 1
        with self.assertRaisesRegex(ValueError, 'protected function changed'):
            self.compare(b,c,calls)

    def test_changed_call_opcode_and_wrong_targets_rejected(self):
        for which, offset in (('b',0),('c',0),('b',1),('c',1)):
            b,c,calls = self.fixture()
            (b if which == 'b' else c)[offset] ^= 1
            with self.assertRaises(ValueError):
                self.compare(b,c,calls)

    def test_missing_duplicate_and_extra_calls_rejected(self):
        b,c,calls = self.fixture()
        for changed in (calls[:3],calls+calls[:1],calls[:2]+calls[1:2]+calls[3:]):
            with self.assertRaises(ValueError):
                self.compare(b,c,changed)


if __name__ == '__main__':
    unittest.main()
