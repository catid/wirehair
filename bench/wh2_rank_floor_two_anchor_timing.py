#!/usr/bin/env python3
"""Fail-closed shim for the retired standalone WH2 timing runner.

Current timing evidence is produced only through the prepared preferred-attempt
controller, which freezes executable identities, bounds process groups, and
validates host isolation.  Git history retains the obsolete ABBA runner; no
part of it remains callable from this compatibility entry point.
"""

from __future__ import annotations

import sys


STANDALONE_RETIREMENT_MESSAGE = (
    "the standalone WH2 rank-floor timing runner is retired and cannot "
    "produce supported evidence; use wh2_preferred_attempt_search.py "
    "run-timing --result-dir <prepared-result-dir>"
)


def main() -> int:
    """Reject the retired runner before inspecting arguments or inputs."""
    print(STANDALONE_RETIREMENT_MESSAGE, file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
