#!/usr/bin/env python3
"""Additive one-shot controller for the .62 GF256 structural screen.

Reuse the reviewed .61 lifecycle without modifying its frozen implementation
or namespace. Importing this wrapper never imports the mathematical worker.
The .62 scientific contract is recorded in wirehair-sxvz.16.1.20.62.
"""

import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "_thue_morse_frozen_controller", HERE / "Wh2NoncommutingRadixRunR0.py")
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)
C.PROTOCOL = "wirehair.wh2.thue-morse-r0"
C.OUTPUT = Path("/var/tmp/wh2-thue-morse-r0")
C.SOURCES = (
    "bench/Wh2ThueMorseR0.py",
    "bench/Wh2ThueMorseRunR0.py",
    "bench/test_Wh2ThueMorseR0.py",
    "bench/test_Wh2ThueMorseRunR0.py",
    "bench/Wh2NoncommutingRadixR0.py",
    "bench/Wh2NoncommutingRadixRunR0.py",
    "bench/test_Wh2NoncommutingRadixR0.py",
    "bench/test_Wh2NoncommutingRadixRunR0.py",
)
_strict_json = C.strict_json


def strict_json(data):
    value = _strict_json(data)
    if (isinstance(value, dict) and value.get("protocol") == C.PROTOCOL and
            value.get("outcome") == "EXHAUSTED"):
        # The existence proof rules exhaustion out. The inherited run() keeps
        # the raw bytes and publishes INVALID when this parser raises.
        raise ValueError("EXHAUSTED is not a Thue-Morse scientific verdict")
    return value


C.strict_json = strict_json


if __name__ == "__main__":
    try:
        sys.exit(C.main())
    except Exception as error:
        print(type(error).__name__ + ": " + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
