#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON=${PYTHON:-python3}

if [[ "$#" -eq 0 ]]; then
    inputs=("-")
else
    inputs=()
    for input in "$@"; do
        if [[ ! -f "$input" ]]; then
            echo "missing dense-count aggregate: $input" >&2
            exit 1
        fi
        inputs+=("$input")
    done
fi

tmp=$(mktemp)
candidate=$(mktemp)
trap 'rm -f "$tmp" "$candidate"' EXIT

"$PYTHON" "$ROOT_DIR/tables/dense_count_validate.py" aggregate \
    "${inputs[@]}" > "$tmp"

awk -F '\t' '
    NF >= 6 {
        ++rows
        n = $1 + 0
        obs = $2 + 0
        min_dense = $3 + 0
        max_dense_value = $4 + 0
        fail = $5 + 0.0
        source = $6

        observations[n] += obs
        if (!(n in min_generated) || min_dense < min_generated[n]) {
            min_generated[n] = min_dense
        }
        if (!(n in dense_count) || max_dense_value > dense_count[n]) {
            dense_count[n] = max_dense_value
            failures[n] = fail
            sources[n] = source
        }
        else if (max_dense_value == dense_count[n]) {
            if (fail < failures[n]) {
                failures[n] = fail
            }
            sources[n] = sources[n] "," source
        }
    }
    END {
        if (rows == 0) {
            print "no dense-count aggregate rows found" > "/dev/stderr"
            exit 1
        }
        for (n in dense_count) {
            printf "%u\t%u\t%.9g\t%u\t%u\t%s\n",
                n, dense_count[n], failures[n],
                observations[n], min_generated[n], sources[n]
        }
    }
' "$tmp" | LC_ALL=C sort -n -k1,1 > "$candidate"

printf "N\tDenseCount\tSelectedFailureRate\tObservations\tMinGeneratedDense\tSourcesAtMax\n"
cat "$candidate"
