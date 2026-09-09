"""Retained K2 native fixture validation, not a new recovery campaign."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

SPEC=importlib.util.spec_from_file_location('_k2_native_tests',Path(__file__).with_name('Wh2K2NativeDataR0.py'))
M=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(M)


class DataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=M.load_report();cls.lookup=M.build_lookup()

    def test_projection_counts_and_deficiency(self):
        data=M.extract(self.report)
        self.assertEqual({k:len(v) for k,v in data.items()},dict(traces=6216,history=56,windows=30,pairs=1539,rows=5274))
        self.assertEqual(sum(p['rank']==1 for p in data['pairs']),5)
        self.assertEqual(sum(p['tail']!=p['B'] for p in data['history']),1)
        self.assertEqual(len({(p['B'],p['tail'],tuple(p['ids'])) for p in data['history']}),16)
        self.assertEqual(len(self.lookup),7168);self.assertEqual(M.C.sha(self.lookup),M.LOOKUP_SHA)
        self.assertEqual(self.lookup[:4],bytes((1,0,0,1)))

    def test_every_emitted_value(self):
        header=M.render(M.extract(self.report),self.lookup).decode('ascii')
        def values(name):
            body=header.split(name,1)[1].split(' = {\n',1)[1].split('};',1)[0]
            return [json.loads(line.rstrip(',').replace('{','[').replace('}',']').replace('u','')) for line in body.splitlines()]
        self.assertEqual(values('Trace kTraces[]'),[[r['B'],r['ids'],r['ranks']] for r in self.report['fresh']+self.report['hard']])
        self.assertEqual(values('Origin kHistory[]'),[[r['b'],r['tail'],len(r['ids']),r['ids']] for r in self.report['inputs']['origins']])
        self.assertEqual(values('kWindows[][6]'),[r['ids'] for r in self.report['seams']])
        pairs=self.report['legacy']+[p for s in self.report['strides'] for p in s['pairs']]
        self.assertEqual(values('Pair kPairs[]'),[[r['ids'],r['rank']] for r in pairs])
        self.assertEqual(values('Row kRows[]'),[[r['id'],r['row']] for r in self.report['evidence']['unique_rows']])
        body=header.split('kLookup[7168] = {\n',1)[1].split('};',1)[0]
        self.assertEqual(bytes(int(v) for v in body.replace('\n','').rstrip(',').split(',')),self.lookup)
        self.assertIn(M.RAW_SHA,header);self.assertIn(M.LOOKUP_SHA,header)

    def test_no_old_bundle_receipt_or_candidate_import(self):
        with mock.patch.object(M.R,'load_report',side_effect=AssertionError('old bundle')), \
             mock.patch.object(M.C,'current_receipt',side_effect=AssertionError('old HEAD')):
            self.assertEqual(M.extract(M.load_report()),M.extract(self.report))
        self.assertNotIn('Wh2K2ThueMorseR0.py',Path(M.__file__).read_text())

    def test_each_tampered_member_rejected(self):
        original=M.C.read_regular
        for member in ('COMPLETE.json',)+M.R.MEMBERS:
            def changed(path,*args,**kwargs):
                raw=original(path,*args,**kwargs);return raw+b' ' if path.name==member else raw
            with self.subTest(member=member),mock.patch.object(M.C,'read_regular',side_effect=changed):
                with self.assertRaisesRegex(ValueError,'identity'):M.load_report()

    def test_invalid_projection(self):
        mutations=[lambda r:r['fresh'].pop(),lambda r:r['hard'].pop(),lambda r:r['history'].pop(),
            lambda r:r['inputs']['origins'].pop(),lambda r:r['fresh'][0]['ranks'].__setitem__(0,1),
            lambda r:r['fresh'][0].__setitem__('B',True),lambda r:r['hard'][0]['ranks'].__setitem__(0,True),
            lambda r:r['hard'][0]['ids'].__setitem__(0,r['hard'][0]['ids'][1]),
            lambda r:r['inputs']['origins'][-2].__setitem__('tail',1280),
            lambda r:r['inputs']['prefixes'][0]['original_widths'].pop(),
            lambda r:r['history'][0].__setitem__('rank',1),lambda r:r['pair'][0][0].__setitem__(1,4),
            lambda r:r['seams'].pop(),lambda r:r['seams'][0].__setitem__('deficient',[[0,1]]),
            lambda r:r['legacy'].pop(),lambda r:r['strides'].pop(),lambda r:r['strides'][0]['pairs'].pop(),
            lambda r:r['strides'][0]['pairs'][0].__setitem__('rank',True),
            lambda r:r['strides'][0]['pairs'][0].__setitem__('determinant',0),
            lambda r:r['strides'][0].__setitem__('failures',0),
            lambda r:r['evidence']['unique_rows'].pop(),lambda r:r['evidence']['unique_rows'][0]['row'].__setitem__(0,False)]
        for i,mutate in enumerate(mutations):
            r=copy.deepcopy(self.report);mutate(r)
            with self.subTest(i=i),self.assertRaises(ValueError):M.extract(r)

    def test_immutable_output_and_lookup(self):
        for lookup in (self.lookup[:-1],bytes([self.lookup[0]^1])+self.lookup[1:]):
            with self.assertRaisesRegex(ValueError,'lookup bytes/hash'):M.render(M.extract(self.report),lookup)
        with tempfile.TemporaryDirectory(prefix='wh2-k2-data-neutral-') as directory:
            path=Path(directory)/'fixture.inc';M.R.write_header(path,b'neutral');M.R.write_header(path,b'neutral')
            with self.assertRaisesRegex(ValueError,'differs'):M.R.write_header(path,b'changed')
            link=Path(directory)/'link.inc';link.symlink_to(path)
            with self.assertRaisesRegex(ValueError,'regular file'):M.R.write_header(link,b'neutral')
            self.assertEqual(path.read_bytes(),b'neutral')


if __name__=='__main__':unittest.main()
