#!/usr/bin/env python3
"""Artifact-policy successor, preserving the R0 experiment and scientific rules.

This private adapter changes artifact identity, not the measured workload or its
decision. A real run still needs a separately frozen, reviewed artifact/argv.
"""
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("_stack_crossover_r1_frozen_r0",
    ROOT / "bench/Wh2PublicStackCrossoverR0.py")
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)

c.CAMPAIGN = "wh2-public-stack-crossover-r1"
c.PREFIX = "wirehair.wh2.public-stack-crossover-r1"
c.OUTPUT = Path("/var/tmp/" + c.CAMPAIGN)
c.SOURCES = ("bench/Wh2PublicStackCrossoverR1.cpp",
    "bench/Wh2InstalledCodeArtifact.cpp") + c.SOURCES[1:]
c.INPUTS = tuple(dict.fromkeys(c.INPUTS + c.SOURCES + (
    "bench/Wh2PublicStackCrossoverR1.py", "bench/test_Wh2PublicStackCrossoverR1.py",
    "bench/Wh2InstalledCodeArtifact.h", "bench/Wh2InstalledCodeArtifactTest.cpp")))
c.HELPER_HASHES.update({
    "bench/Wh2PublicStackCrossoverR0.cpp": "8bd8dcb19d7c8106bb12bee591982b22d9e7cba5de303aaccb95d5528b2b1144",
    "bench/Wh2PublicStackCrossoverR0.py": "c456ce01af390cd65ea5db498dfce79ab9e4d8bba3f3d70150b76d6ef6c05178",
    "bench/test_Wh2PublicStackCrossoverR0.py": "d672985dfe65de99a401c6cde90b7f297b7705e35668166d3693f3c27ae5743d",
})
c.h.CAMPAIGN, c.h.FIXED_OUTPUT_DIR, c.h.R1_TRACKED_INPUTS = c.CAMPAIGN, c.OUTPUT, c.INPUTS
c.h.R0_ROOTS += (Path("/var/tmp/wh2-public-stack-crossover-r0"),)

# Verify the real installed library using the exact worker before any namespace
# claim and again at postflight. This probe has no dlopen, codec calls or clocks.
# Bind the override in the private module: its build/run functions use its globals.
_frozen_check_build = c.check_build


def check_build(receipt, live=True):
    _frozen_check_build(receipt, live=live)
    if live:
        c.equal(c.command([receipt["worker_path"], "--artifact-preflight"]),
            "stack-crossover-r1 installed libc verified (no codec or clocks)\n",
            "installed-code preflight")


c.check_build = check_build

if __name__ == "__main__":
    try:
        sys.exit(c.main())
    except (c.h.ValidationError, OSError, ValueError, TypeError, KeyError, IndexError) as error:
        print("Stack crossover r1: " + str(error), file=sys.stderr)
        sys.exit(1)
