#!/usr/bin/env python3
"""R1: same preserved-path gate, with positive worker claim-path qualification.

R0 is spent INVALID before any codec work because its worker read the original
admission claim path. No R0 samples or namespace are reused. All WORK, cases,
timing schedules, controls, thresholds, bounds and libraries remain unchanged.
"""
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location('small_isolation_r0_definition',
    Path(__file__).with_name('Wh2SmallIsolationPreservedCostR0.py'))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
SETTINGS = M.SETTINGS._replace(
    protocol='wirehair.wh2.small-isolation-preserved-cost-r1',
    output=Path('/var/tmp/wh2-small-isolation-preserved-cost-r1'),
    sources=M.NEW+('bench/Wh2SmallIsolationPreservedCostR1.py',
                   'bench/test_Wh2SmallIsolationPreservedCostR1.py'))

if __name__ == '__main__': M.main(SETTINGS)
