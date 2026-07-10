#!/usr/bin/env python3
"""Reject disproven heavy-matrix guarantees from source documentation."""

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKIP_PARTS = {".git", ".beads", ".claude", ".codex", "__pycache__"}
SUFFIXES = {".c", ".cc", ".cpp", ".h", ".hpp", ".md", ".txt"}
PROHIBITED = (
    r"any\s+disturbance\s+from\s+the\s+binary\s+rows\s+above\s+does\s+not\s+affect",
    r"100%\s+invertibility\s+of\s+the\s+matrix",
    r"heavy.{0,120}guaranteed\s+pivots",
    r"heavy\s+rows?.{0,120}guaranteed\s+invertib",
)


def source_files():
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if (not path.is_file() or path.suffix.lower() not in SUFFIXES or
                any(part in SKIP_PARTS or part == "build" or
                    part.startswith("build-") for part in relative.parts)):
            continue
        yield path


def main():
    failures = []
    for path in source_files():
        text = " ".join(path.read_text(encoding="utf-8", errors="replace").split()).lower()
        for pattern in PROHIBITED:
            if re.search(pattern, text):
                failures.append("%s: matches %s" % (path.relative_to(ROOT), pattern))

    authority = (ROOT / "tables" / "HEAVY_MATRIX.md").read_text(encoding="utf-8").lower()
    for phrase in ("wire format", "empirical", "1/256", "not a reliability guarantee",
                   "additional recovery rows"):
        if phrase not in authority:
            failures.append("tables/HEAVY_MATRIX.md: missing %r" % phrase)
    for name in ("GenerateDenseCount.cpp", "GeneratePeelSeeds.cpp",
                 "GenerateMostDenseSeeds.cpp", "GenerateSmallDenseSeeds.cpp"):
        text = (ROOT / "tables" / name).read_text(encoding="utf-8")
        if "HEAVY_MATRIX.md" not in text:
            failures.append("tables/%s: missing authoritative-doc reference" % name)

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("heavy-matrix documentation lint passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
