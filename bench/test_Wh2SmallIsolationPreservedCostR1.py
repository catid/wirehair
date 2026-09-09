#!/usr/bin/env python3
"""Configuration/binding tests; never launch a scientific namespace."""
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location('small_isolation_r1_tested',
    Path(__file__).with_name('Wh2SmallIsolationPreservedCostR1.py'))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class Tests(unittest.TestCase):
    def test_separate_namespace_without_changing_work(self):
        previous=M.M.SETTINGS; current=M.SETTINGS
        self.assertNotEqual(previous.output,current.output)
        self.assertNotEqual(previous.protocol,current.protocol)
        self.assertEqual(previous.libraries,current.libraries)
        self.assertIs(previous.provenance,current.provenance)
        self.assertEqual(current.sources[:-2],previous.sources)

    def test_worker_claim_binding_is_not_the_spent_path(self):
        record=M.M.R.claim_binding(M.SETTINGS)
        self.assertEqual(record,dict(protocol='wirehair.wh2.small-isolation-preserved-cost-r1',
            claim_path='/var/tmp/wh2-small-isolation-preserved-cost-r1/CLAIM.json',
            positive_returncode=0,negative_returncode=1))
        for cfg in (M.M.SETTINGS,M.M.R.configuration(None)):
            self.assertNotEqual(M.M.R.claim_binding(cfg)['claim_path'],record['claim_path'])


if __name__=='__main__': unittest.main()
