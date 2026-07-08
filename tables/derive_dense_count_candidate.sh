#!/usr/bin/env bash
set -euo pipefail

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

read_aggregate() {
    local input="$1"
    local source="$input"
    if [[ "$input" == "-" ]]; then
        source="<stdin>"
    fi

    awk -F '\t' -v source="$source" '
        function is_uint(value) {
            return value ~ /^[0-9]+$/
        }

        function is_number(value) {
            return value ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
        }

        FNR == 1 && $1 == "N" { next }
        NF == 0 { next }
        NF >= 5 && is_uint($1) && is_uint($2) && is_uint($3) &&
            is_uint($4) && is_number($5) {
            if (($3 + 0) > ($4 + 0)) {
                printf "malformed dense-count aggregate row: %s:%u: min > max: %s\n",
                    source, FNR, $0 > "/dev/stderr"
                exit 1
            }
            row_source = (NF >= 6 && $6 != "") ? $6 : source
            printf "%s\t%s\t%s\t%s\t%s\t%s\n",
                $1, $2, $3, $4, $5, row_source
            next
        }
        {
            printf "malformed dense-count aggregate row: %s:%u: %s\n",
                source, FNR, $0 > "/dev/stderr"
            exit 1
        }
    ' "$input"
}

for input in "${inputs[@]}"; do
    read_aggregate "$input"
done > "$tmp"

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

printf "N\tDenseCount\tLowestFailures\tObservations\tMinGeneratedDense\tSourcesAtMax\n"
cat "$candidate"
