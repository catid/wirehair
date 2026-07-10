#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RESULT_DIR="${1:-$ROOT/experiments/precode/results}"
PYTHON=${PYTHON:-python3}
VALIDATOR="$ROOT/experiments/precode/validate_results.py"

if [ ! -d "$RESULT_DIR" ]; then
  echo "missing results directory: $RESULT_DIR" >&2
  exit 1
fi

fail=0

mapfile -d '' csv_files < <(find "$RESULT_DIR" -type f -name '*.csv' -print0)
if ((${#csv_files[@]} > 0)); then
  if ! "$PYTHON" "$VALIDATOR" "${csv_files[@]}"; then
    fail=1
  fi
else
  echo "no CSV result files: $RESULT_DIR" >&2
  fail=1
fi

while IFS= read -r err; do
  if [ -s "$err" ]; then
    echo "non-empty error log: $err" >&2
    fail=1
  fi
done < <(find "$RESULT_DIR" -type f -name '*.err' | sort)

while IFS= read -r tmp; do
  echo "stale temporary result file: $tmp" >&2
  fail=1
done < <(find "$RESULT_DIR" -type f \( -name '*.tmp.*' -o -name '*.err.tmp.*' \) | sort)

if [ "$fail" -ne 0 ]; then
  exit 1
fi

echo "precode result integrity OK: $RESULT_DIR"
