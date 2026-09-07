#!/usr/bin/env python3
"""Additive .82 launcher, reusing the unchanged reviewed .61 capture lifecycle.

Receipt mode only authenticates committed sources and historical bytes. It
never imports the new mathematical worker or evaluates any candidate.
"""
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C = sibling("_k3_tm_controller", "Wh2NoncommutingRadixRunR0.py")
H = sibling("_k3_tm_receipt_history", "Wh2ThueMorseRecoveryHistoryR0.py")
C.PROTOCOL = "wirehair.wh2.k3-thue-morse-r0"
C.OUTPUT = Path("/var/tmp/wh2-k3-thue-morse-r0")
C.SOURCES = (
    "bench/Wh2K3ThueMorseR0.py",
    "bench/Wh2K3ThueMorseRunR0.py",
    "bench/test_Wh2K3ThueMorseR0.py",
    "bench/Wh2NoncommutingRadixR0.py",
    "bench/Wh2NoncommutingRadixRunR0.py",
    "bench/test_Wh2NoncommutingRadixRunR0.py",
    "bench/Wh2ThueMorseRecoveryHistoryR0.py",
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h",
    "bench/wh2_benchmark_contract_v4.json",
)
_current_receipt = C.current_receipt


def current_receipt(deadline=None):
    receipt = _current_receipt(deadline)
    _, provenance = H.read_bundle(H.WIDTH_LOCAL, "COMPLETE", H.WIDTH_MANIFEST,
                                  H.WIDTH_FILES, False, [0], deadline)
    receipt["history"] = provenance
    return receipt


C.current_receipt = current_receipt


if __name__ == "__main__":
    try:
        sys.exit(C.main())
    except Exception as error:
        print(type(error).__name__ + ": " + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
