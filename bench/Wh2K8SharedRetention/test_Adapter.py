"""Synthetic/source-only gate-A tests; never run codecs or spent controllers."""
import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import Adapter as K


class AdapterTest(unittest.TestCase):
    @staticmethod
    def records(ratio=1):
        return [dict(coordinate=c, observation=dict(clocks=[0,0,0,
            round(1000000*(ratio if c[5] == 2 and K.R.SIDES[c[6]] ^ c[2] else 1)),0,0]))
            for c in K.R.roster(K.CASES)]

    def test_source_derivation_and_unchanged_work(self):
        common, worker = K.worker_sources()
        original = (K.ROOT/'bench/Wh2AdmissionRegressionCostR0.cpp').read_bytes()
        def work(raw): return raw[raw.index(b'NOINLINE void RunWork('):raw.index(b'void Prepare(')]
        self.assertEqual(work(common), work(original))
        self.assertIn(b'case_count=38,max_batch=128,callbacks=98496', common)
        self.assertIn(b'max_batch*(22*1280+128)', common)
        self.assertIn(b'#include "Wh2K8RetentionCommon.h"', worker)
        self.assertNotIn(b'slot_count', common)
        for raw in (common, worker):
            self.assertIn(b'cpu={210,210}', raw)
            self.assertNotIn(b'UINT64_C(180000000000)', raw)
            self.assertIn(b'UINT64_C(150000000000)', raw)

    def test_exact_roster_without_mutating_historical_defaults(self):
        self.assertEqual(K.CASES, K.F.CASES)
        self.assertEqual(len(set(K.CASES)), 38)
        self.assertEqual(K.R.callback_count(K.CASES), 98496)
        self.assertEqual(K.R.callback_count(), 51840)
        self.assertEqual(K.D.SETTINGS.cases, K.R.CASES)
        self.assertNotEqual(K.PROTOCOL, K.F.PROTOCOL)
        self.assertNotEqual(K.OUTPUT, K.F.OUTPUT)
        self.assertNotIn(K.F.CANDIDATE_SHA, [sha for _,sha in K.LIBRARIES])

    def test_retention_needs_no_incidental_certified_speedup(self):
        passed = dict(outcome='PASS')
        for outcome in ('PASS', 'CONTROL_FAIL', 'REGRESSION', 'INCONCLUSIVE'):
            result = K.combine([passed, dict(outcome=outcome)])
            self.assertEqual(result['outcome'], outcome)
            self.assertEqual(result['current_path_retention_qualified'], outcome == 'PASS')
            for name in ('pre_admission_restoration_qualified', 'ordinary_K8_WH1_speed_qualified',
                         'historical_K3_workload_retention_qualified', 'production_promotion_claimed'):
                self.assertFalse(result[name])
        self.assertNotIn('certified_four_cell_improvements', K.combine([passed, passed]))

    def test_both_load_orders_are_required(self):
        with self.assertRaises(ValueError): K.combine([dict(outcome='PASS')])

    def test_every_outcome_pair_preserves_failure_precedence(self):
        outcomes = ('CONTROL_FAIL', 'REGRESSION', 'INCONCLUSIVE', 'PASS')
        for first in outcomes:
            for second in outcomes:
                result = K.combine([dict(outcome=first), dict(outcome=second)])
                expected = min((first, second), key=outcomes.index)
                self.assertEqual(result['outcome'], expected)
                self.assertEqual(result['current_path_retention_qualified'], expected == 'PASS')
        with self.assertRaises(ValueError):
            K.combine([dict(outcome='PASS'), dict(outcome='UNKNOWN')])

    def test_configuration_routes_only_new_explicit_adapter(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k8-config-test.') as directory:
            adapter = Path(directory)/'adapter'
            K.prepare(adapter)
            cfg = K.settings(adapter)
            self.assertEqual((cfg.protocol, cfg.output, cfg.libraries, cfg.cases),
                             (K.PROTOCOL, K.OUTPUT, K.LIBRARIES, K.CASES))
            self.assertIs(cfg.header_checker, K.verify_header)
            self.assertIs(cfg.result_combiner, K.combine)
            self.assertIs(cfg.qualify, K.qualify)
            self.assertEqual(cfg.worker_source, str(adapter/K.WORKER))
            self.assertIn(cfg.worker_source, cfg.sources)
            self.assertEqual(K.R.claim_binding(cfg)['claim_path'], str(K.OUTPUT/'CLAIM.json'))
            with patch.object(K, 'provenance', return_value=('new-proof', set())) as proof:
                self.assertEqual(cfg.provenance(), ('new-proof', set()))
                proof.assert_called_once_with(adapter.resolve(), None)
                proof.reset_mock()
                cfg.provenance(Path(directory)/'new-build')
                proof.assert_called_once_with(adapter.resolve(), Path(directory)/'new-build')

    def test_complete_new_roster_and_every_phase(self):
        roster = list(K.R.roster(K.CASES))
        seen, phases = {}, {}
        for index, c in enumerate(roster):
            i,r,o,w,m,pair,p,arm,delay = c
            self.assertEqual(i, index)
            self.assertEqual(arm, K.R.PAIRS[pair][K.R.SIDES[p] ^ o])
            key = w,m,pair,o
            if p == 0: seen.setdefault(key, []).append(r)
            if p >= 2 and p % 2 == 0:
                self.assertEqual(delay, roster[index+1][-1])
                self.assertNotEqual(K.R.SIDES[p], K.R.SIDES[p+1])
                phases.setdefault(key, []).append((delay*96//1000000)//2)
        self.assertEqual(len(roster), 98496)
        self.assertEqual(len(seen), 456)
        self.assertEqual(sum(key[2] < 2 for key in seen), 304)
        for key in seen:
            self.assertEqual(sorted(seen[key]), list(range(12)))
            self.assertEqual(sorted(phases[key]), sorted(list(range(48))*2))

    def test_new_roster_retention_statistics_and_control_dominance(self):
        for ratio, outcome in ((1, 'PASS'), (.999, 'PASS'), (1.001, 'REGRESSION')):
            result = K.R.statistics(self.records(ratio), K.CASES)
            self.assertEqual(result['outcome'], outcome)
            self.assertEqual(len(result['statistics']), 456)
            self.assertEqual(len(result['resolved_regressions']), 152 if ratio > 1 else 0)
            self.assertEqual(K.combine([result, result])['current_path_retention_qualified'], outcome == 'PASS')
        rows = self.records(.8)
        for row in rows:
            c = row['coordinate']
            if c[3:6] == [37,1,1] and K.R.SIDES[c[6]] ^ c[2]:
                row['observation']['clocks'][3] = 1100000
        result = K.R.statistics(rows, K.CASES)
        self.assertEqual(result['outcome'], 'CONTROL_FAIL')
        self.assertEqual(len(result['failed_controls']), 2)
        with self.assertRaises(ValueError): K.R.statistics(rows[:-1], K.CASES)
        with self.assertRaises(ValueError): K.R.statistics(rows)
        rows = self.records()
        for row in rows:
            c = row['coordinate']
            if c[3:6] == [37,1,2] and K.R.SIDES[c[6]] ^ c[2]:
                row['observation']['clocks'][3] = 1001000
        result = K.R.statistics(rows, K.CASES)
        self.assertEqual(result['resolved_regressions'], [[37,1,2,0], [37,1,2,1]])
        self.assertFalse(K.combine([result, result])['current_path_retention_qualified'])

    def test_full_raw_uses_actual_endpoint_and_rejects_late_corruption(self):
        zero = [0]*4
        previous = dict(clocks=[1,0,2,3,0,4], before=zero, after=zero)
        # Deliberately different arm endpoints catch accidental K-only or
        # shared-endpoint accounting. This is synthetic data, not a codec claim.
        header = dict(prelude=previous, fixtures=[dict(arms=[dict(steps=c[1]),
            dict(steps=c[1]+1)]) for c in K.CASES])
        raw = [header]
        addresses = {n:list(range(1,n+1))+[0]*(128-n) for n in (4,128)}
        for c in K.R.roster(K.CASES):
            ready = previous['clocks'][5]+1; cpu = previous['clocks'][4]
            target = ready+c[-1]; start = target+3; case = K.CASES[c[3]]
            cycles, metric = K.R.batch(case), c[4]
            steps = header['fixtures'][c[3]]['arms'][c[7]]['steps']
            observed = dict(clocks=[start,cpu+2,start+1,start+1000001,cpu+1000002,start+1000002],
                            before=zero, after=zero)
            raw.append(dict(type='record', coordinate=c, ready=ready, target=target,
                wait=[ready+1,cpu,target+2,cpu+1], observation=observed,
                counts=[0 if metric else cycles,0 if metric else cycles*(case[1]+14),
                        cycles if metric else 0,cycles*steps if metric else 0,
                        cycles if metric else 0,cycles],
                addresses=addresses[cycles], address_count=cycles, complete=True, checked=True))
            previous = observed
        count = K.R.callback_count(K.CASES)
        raw.append(dict(type='footer', complete=True, records=count, work_ns=count*1000000))
        with patch.object(K, 'verify_header') as checker:
            def verify():
                return K.R.verify_rows(raw, 'claim', 1, [], K.PROTOCOL, K.CASES, checker)
            self.assertEqual(verify()['outcome'], 'PASS')
            checker.assert_called_once_with(header, 1, 'claim', [], K.PROTOCOL)
            index = next(i for i in range(1,len(raw)-1)
                         if raw[i]['coordinate'][3:5] == [37,1] and raw[i]['coordinate'][7] == 1)
            saved = raw[index]
            bad_counts = list(saved['counts']); bad_counts[3] -= 128
            for field, value in (('counts', bad_counts), ('complete', False), ('checked', False),
                                 ('target', 0), ('address_count', 4)):
                raw[index] = dict(saved, **{field:value})
                with self.assertRaises(ValueError): verify()
            raw[index] = saved
            raw[-1]['records'] -= 1
            with self.assertRaisesRegex(ValueError, 'terminal complete footer'): verify()

    def test_adapter_proof_rejects_source_or_helper_drift(self):
        with tempfile.TemporaryDirectory(prefix='wh2-k8-observer-test.') as directory:
            path = Path(directory)/'adapter'
            K.prepare(path)
            self.assertEqual(len(K.verify_adapter(path)), 3)
            with patch.object(K.S, 'worker_sources', return_value=(b'changed', b'changed')):
                with self.assertRaises(ValueError): K.verify_adapter(path)
            with patch.dict(K.HELPERS, {'Wh2AdmissionRegressionCostR0.py':'0'*64}):
                with self.assertRaises(ValueError): K.verify_adapter(path)

    def test_explicit_metadata_rejects_five_slot_or_wrong_library(self):
        meta = [dict(path=str(path), sha256=digest, context_bytes=141328,
                     slots=[dict(name=name) for name in K.SLOTS],
                     exports=[dict(name='symbol'+str(i)) for i in range(53)],
                     runtime_slots=[{}]*37) for path, digest in K.LIBRARIES]
        K.check_metadata(meta)
        for field in ('path', 'sha256', 'context_bytes', 'slots', 'exports', 'runtime_slots'):
            bad = copy.deepcopy(meta)
            bad[1][field] = bad[1][field][:-1] if isinstance(bad[1][field], (list,str)) else 137232
            with self.assertRaises(ValueError): K.check_metadata(bad)
        for bad in (meta[::-1], [meta[0], meta[0]]):
            with self.assertRaises(ValueError): K.check_metadata(bad)
        bad = copy.deepcopy(meta)
        bad[1]['slots'][0], bad[1]['slots'][1] = bad[1]['slots'][1], bad[1]['slots'][0]
        with self.assertRaises(ValueError): K.check_metadata(bad)
        bad = copy.deepcopy(meta)
        bad[1]['path'] = str(K.S.PREPARED/'libwirehair.so.2.0.0')
        bad[1]['sha256'] = K.F.CANDIDATE_SHA
        with self.assertRaises(ValueError): K.check_metadata(bad)

    def test_crossed_bindings_and_overlapping_contexts_rejected(self):
        exports = [dict(name=name, offset=64*(i+1)) for i,name in enumerate(K.SLOTS)]
        meta = [dict(context=4096, context_bytes=141328, exports=exports,
                     slots=[dict(name=name) for name in K.SLOTS]) for _ in range(2)]
        bindings = []
        for i,lib in enumerate(meta):
            base = (i+1)*0x1000000
            targets = [base+s['offset'] for s in exports]
            bindings.append(dict(base=base, context=base+lib['context'], exports=targets,
                slots=list(targets), providers=[10,20,30,40,50], runtime_targets=list(range(1,38))))
        K.R.verify_bindings(bindings, meta)
        for field in ('slots', 'exports', 'runtime_targets', 'providers'):
            bad = copy.deepcopy(bindings)
            bad[1][field][0] = bindings[0][field][0] if field in ('slots','exports') else 999
            with self.assertRaises(ValueError): K.R.verify_bindings(bad, meta)
        with self.assertRaisesRegex(ValueError, 'distinct GF storage'):
            K.R.verify_bindings([bindings[0], bindings[0]], meta)


if __name__ == '__main__': unittest.main()
