from pathlib import Path
import tempfile
import unittest
from unittest import mock
import Generate as G


class Test(unittest.TestCase):
    def test_generated_contract(self):
        source = G.generate()
        for required in ('Batch = 64, Reps = 12','route < 7;',
                         'ks[7] = {3,3,5,8,3,16,3}',
                         'fixtures.size() == 70','Reps*210*5*2*18',
                         'exact standalone WHK3 descriptor',
                         'api.standalone = matched.shape.route == 6 && !api.wh1;',
                         'Work(api,matched,metric,nullptr,output.data()+16,produced);',
                         'WirehairSmall_BorrowedImmutable : WirehairSmall_Independent',
                         'lib->sfree(static_cast<WirehairSmallCodec>(h))',
                         'CPU_COUNT(&mask) == 1','sanitized observer is neutral only',
                         '240-second worker cap includes publication'):
            self.assertIn(required,source)
        self.assertNotIn('Batch = 32',source)
        self.assertEqual(source.count('210'),4)

    def test_exact_sites(self):
        self.assertEqual(G.replace_once('abc','b','d'),'adc')
        for source in ('none','twice twice'):
            with self.assertRaises(ValueError): G.replace_once(source,'twice','x')

    def test_shared_source_hash(self):
        with tempfile.TemporaryDirectory() as temporary, mock.patch.object(G,'ROOT',Path(temporary)):
            path = Path(temporary)/'source'
            path.write_bytes(b'tampered')
            with self.assertRaisesRegex(ValueError,'review shared source changes'):
                G.pinned('source','0'*64)


if __name__ == '__main__':
    unittest.main()
