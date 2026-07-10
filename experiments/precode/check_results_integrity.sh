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
validated_files=0
validated_rows=0

mapfile -d '' csv_files < <(find "$RESULT_DIR" -type f -name '*.csv' -print0)
if ((${#csv_files[@]} > 0)); then
  # This directory aggregates independent campaigns whose control rows may
  # intentionally overlap.  Validate each artifact here; queue promotion
  # separately validates all shards belonging to the campaign being extended.
  for csv_file in "${csv_files[@]}"; do
    if validation_output=$("$PYTHON" "$VALIDATOR" "$csv_file"); then
      if [[ "$validation_output" =~ ,\ ([0-9]+)\ row\(s\)$ ]]; then
        ((validated_files += 1))
        ((validated_rows += BASH_REMATCH[1]))
      else
        echo "unexpected validator output for $csv_file: $validation_output" >&2
        fail=1
      fi
    else
      fail=1
    fi
  done
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

echo "validated $validated_files precode result file(s), $validated_rows row(s)"
echo "precode result integrity OK: $RESULT_DIR"
