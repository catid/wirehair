"""Synthetic and read-only qualification tests; no scientific worker launches."""
import copy
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('current_preserved_tested',
    Path(__file__).with_name('Wh2CurrentPreservedDeferredCostR0.py'))
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)


def fixture(mode):
    previous=dict(clocks=[1,0,2,3,0,4],before=[0]*4,after=[0]*4)
    header=dict(prelude=previous,fixtures=[dict(arms=[dict(steps=2)]*2)],bindings=[])
    rows=[header]; total=0
    for i,coordinate in enumerate(list(C.R.roster())[54:57]):
        ready=previous['clocks'][5]; cpu=previous['clocks'][4]; start=ready+1
        final=i==2
        clocks=[start,cpu+1,start+10,start+20,cpu+2,start+30]
        if final and mode=='last-clock': clocks[3:]=[0,0,0]
        else: total+=10
        obs=dict(clocks=clocks,before=[0]*4,after=[0]*4)
        rows.append(dict(type='record',coordinate=coordinate,ready=ready,target=ready,
            wait=[ready,cpu,ready,cpu],observation=obs,counts=[0,0,128,256,128,128],
            addresses=[123]*128,address_count=128,
            complete=not(final and mode in ('last-recover','throw-recover')),
            checked=not final or mode=='success'))
        previous=obs
    rows.append(dict(type='footer',complete=mode=='success',records=3,work_ns=total))
    return header,rows


def raw(rows): return b''.join(C.A.canonical(row) for row in rows)


class Tests(unittest.TestCase):
    def test_work_symbol_allows_cleanup_lambda_not_duplicate_work(self):
        name=('(anonymous namespace)::RunWork((anonymous namespace)::Api const&, '
              '(anonymous namespace)::Fixture const&, (anonymous namespace)::Arm const&, '
              'unsigned int, unsigned char*, (anonymous namespace)::Work&)')
        body='000000000041297e t '+name
        helper='0000000000407742 t '+name+'::{lambda()#1}::operator()() const'
        C.R.verify_work_symbol(body+'\n'+helper)
        for invalid in (helper,body+'\n'+body,'',body+' [clone .cold]'):
            with self.assertRaisesRegex(ValueError,'single common WORK'):
                C.R.verify_work_symbol(invalid)

    def test_explicit_adapter_and_unchanged_defaults(self):
        old=C.R.configuration(None)
        self.assertIsNone(old.worker_source); self.assertIsNone(old.qualify)
        self.assertEqual(old.protocol,C.R.PROTOCOL)
        self.assertEqual(old.sources,C.R.NEW)
        self.assertNotEqual(old.output,C.OUTPUT)
        self.assertEqual(C.SETTINGS.libraries[0],old.libraries[0])
        self.assertNotEqual(C.SETTINGS.libraries[1],old.libraries[1])
        self.assertIn(C.SETTINGS.worker_source,C.SETTINGS.sources)
        self.assertIs(C.SETTINGS.qualify,C.qualify)
        self.assertEqual(len(C.R.CASES),20)
        self.assertEqual(len(list(C.R.roster())),51840)

    def test_actual_current_metadata(self):
        metadata=C.R.metadata(C.LIBRARIES)
        for m in metadata:
            self.assertEqual((len(m['exports']),len(m['slots']),len(m['runtime_slots'])),(53,6,37))
            self.assertEqual(m['context_bytes'],0x22810)
        self.assertEqual(metadata[1]['sha256'],C.LIBRARIES[1][1])

    def test_synthetic_publication_all_outcomes(self):
        with patch.object(C.R,'verify_header'):
            for mode in C.PUBLICATIONS:
                header,rows=fixture(mode)
                self.assertEqual(C.verify_publication(raw(rows),mode,0,[],header),
                    dict(mode=mode,load_order=0,records=3,scientific_launch=False))

    def test_every_record_prefix_and_failure_field_checked(self):
        for mode in C.PUBLICATIONS:
            header,rows=fixture(mode)
            with patch.object(C.R,'verify_header'):
                for mutate in (lambda r:r.pop(),lambda r:r[3].__setitem__('checked',not r[3]['checked']),
                    lambda r:r[1]['counts'].__setitem__(5,127),lambda r:r[3]['addresses'].__setitem__(127,0),
                    lambda r:r[3].__setitem__('complete',not r[3]['complete']),
                    lambda r:r[2].__setitem__('target',0),lambda r:r[-1].__setitem__('work_ns',0)):
                    changed=copy.deepcopy(rows); mutate(changed)
                    with self.assertRaises(ValueError): C.verify_publication(raw(changed),mode,0,[],header)

    def test_partial_clock_not_silently_completed(self):
        header,rows=fixture('last-clock'); rows[3]['observation']['clocks'][3]=999
        with patch.object(C.R,'verify_header'):
            with self.assertRaisesRegex(ValueError,'unfinished final clock'):
                C.verify_publication(raw(rows),'last-clock',0,[],header)

    def test_initial_fixture_header_is_not_replaced_by_late_mutation(self):
        header,rows=fixture('last-source'); header=copy.deepcopy(header)
        rows[0]['fixtures'][0]['arms'][1]['steps']=3
        with patch.object(C.R,'verify_header'):
            with self.assertRaisesRegex(ValueError,'immutable initial header'):
                C.verify_publication(raw(rows),'last-source',0,[],header)

    def test_only_binding_bases_may_vary_between_neutral_processes(self):
        header,rows=fixture('success'); header=copy.deepcopy(header); rows[0]['bindings']=[123]
        with patch.object(C.R,'verify_header') as verify:
            C.verify_publication(raw(rows),'success',1,[],header)
            verify.assert_called_once()

    def test_claim_namespace_is_fresh_and_separate(self):
        claim=C.R.claim_binding(C.SETTINGS)
        self.assertEqual(claim['claim_path'],str(C.OUTPUT/'CLAIM.json'))
        self.assertNotEqual(claim,C.R.claim_binding(C.R.configuration(None)))
        self.assertNotEqual(C.SETTINGS.protocol,C.H.PROTOCOL)

    def test_original_work_source_is_included_not_copied(self):
        source=(C.ROOT/C.SETTINGS.worker_source).read_text()
        self.assertIn('#include "Wh2AdmissionRegressionCostR0.cpp"',source)
        self.assertNotIn('void RunWork(',source)
        body=source[source.index('int DeferredWorker('):source.index('struct PublicationReader')]
        self.assertNotIn('RecordJson(',body)
        self.assertNotIn('HeaderJson(',body)
        self.assertEqual(body.count('PublishDeferred('),1)


if __name__=='__main__': unittest.main()
