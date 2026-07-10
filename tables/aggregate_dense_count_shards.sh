#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON=${PYTHON:-python3}

inputs=()
if [[ "$#" -gt 0 ]]; then
    for input in "$@"; do
        if [[ -d "$input" ]]; then
            input="$input/dense_count.out"
        fi
        if [[ ! -f "$input" ]]; then
            echo "missing dense-count shard: $input" >&2
            exit 1
        fi
        inputs+=("$input")
    done
else
    while IFS= read -r input; do
        inputs+=("$input")
    done < <(find "$ROOT_DIR/tables/results" -maxdepth 2 \
        -name dense_count.out -type f 2>/dev/null | sort)
fi

if [[ "${#inputs[@]}" -eq 0 ]]; then
    echo "usage: $0 [dense_count.out-or-shard-dir ...]" >&2
    echo "or place shard directories under tables/results/" >&2
    exit 1
fi

tmp=$(mktemp)
aggregate=$(mktemp)
trap 'rm -f "$tmp" "$aggregate"' EXIT

"$PYTHON" "$ROOT_DIR/tables/dense_count_validate.py" shards \
    "${inputs[@]}" > "$tmp"

awk -F '\t' '
    NF >= 4 {
        ++rows
        n = $1 + 0
        dense = $2 + 0
        fail = $3 + 0.0
        source = $4
        ++obs[n]
        if (!(n in min_dense) || dense < min_dense[n]) {
            min_dense[n] = dense
        }
        if (!(n in max_dense) || dense > max_dense[n]) {
            max_dense[n] = dense
            fail_at_max[n] = fail
            sources_at_max[n] = source
        }
        else if (dense == max_dense[n]) {
            if (fail < fail_at_max[n]) {
                fail_at_max[n] = fail
            }
            sources_at_max[n] = sources_at_max[n] "," source
        }
    }
    END {
        if (rows == 0) {
            print "no dense-count rows found" > "/dev/stderr"
            exit 1
        }
        for (n in obs) {
            printf "%u\t%u\t%u\t%u\t%.9g\t%s\n",
                n, obs[n], min_dense[n], max_dense[n],
                fail_at_max[n], sources_at_max[n]
        }
    }
' "$tmp" | LC_ALL=C sort -n -k1,1 > "$aggregate"

printf "N\tObservations\tMinGeneratedDense\tMaxGeneratedDense\tLowestFailuresAtMax\tSourcesAtMax\n"
cat "$aggregate"
