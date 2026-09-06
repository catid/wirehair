#!/usr/bin/env python3
"""One-shot .63 rank feasibility launcher; no codec or timing claims.

The inherited lifecycle keeps historical controllers and namespaces unchanged.
Receipt creation authenticates existing evidence but never imports the worker,
derives fresh roots, or scores the fixed candidate. The scientific protocol is
preregistered in wirehair-sxvz.16.1.20.63.
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


C = sibling("_thue_recovery_controller", "Wh2NoncommutingRadixRunR0.py")
H = sibling("_thue_recovery_history_receipts", "Wh2ThueMorseRecoveryHistoryR0.py")
C.PROTOCOL = "wirehair.wh2.thue-morse-recovery-r0"
C.OUTPUT = Path("/var/tmp/wh2-thue-morse-recovery-r0")
C.SOURCES = (
    "bench/Wh2ThueMorseRecoveryR0.py",
    "bench/Wh2ThueMorseRecoveryRunR0.py",
    "bench/test_Wh2ThueMorseRecoveryR0.py",
    "bench/Wh2ThueMorseRecoveryHistoryR0.py",
    "bench/test_Wh2ThueMorseRecoveryHistoryR0.py",
    "bench/test_Wh2ThueMorseRecoveryRunR0.py",
    "bench/Wh2ThueMorseR0.py",
    "bench/Wh2NoncommutingRadixR0.py",
    "bench/Wh2NoncommutingRadixRunR0.py",
    "bench/test_Wh2ThueMorseR0.py",
    "bench/test_Wh2NoncommutingRadixR0.py",
    "bench/test_Wh2NoncommutingRadixRunR0.py",
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h",
    "bench/Wh2FrozenTraceTest.cpp",
    "bench/wh2_benchmark_contract_v4.json",
    "bench/wh2_mix2_seed_repair_contract.py",
    "bench/wh2_precodefail_work_screen.py",
)
_current_receipt = C.current_receipt
_strict_json = C.strict_json


def current_receipt(deadline=None):
    receipt = _current_receipt(deadline=deadline)
    receipt["historical_inputs"] = H.input_receipt(deadline=deadline)
    return receipt


def strict_json(data):
    value = _strict_json(data)
    if (isinstance(value, dict) and value.get("protocol") == C.PROTOCOL and
            value.get("outcome") == "EXHAUSTED"):
        raise ValueError("EXHAUSTED is not a fixed-candidate recovery verdict")
    return value


C.current_receipt = current_receipt
C.strict_json = strict_json


if __name__ == "__main__":
    try:
        sys.exit(C.main())
    except Exception as error:
        print(type(error).__name__ + ": " + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
