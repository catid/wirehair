#!/bin/bash
# Validate the legacy seedsearch positional contract before a generator writes output.
set -euo pipefail

if [ "$#" -eq 0 ]; then
    echo "usage: validate_seedsearch.sh FILE..." >&2
    exit 2
fi

for input in "$@"; do
    if [ ! -r "$input" ]; then
        echo "seedsearch input is not readable: $input" >&2
        exit 2
    fi
    awk -v file="$input" '
        function uint(x) { return x ~ /^[0-9]+$/ }
        function sint(x) { return x ~ /^-?[0-9]+$/ }
        function num(x) {
            return x ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/
        }
        function score(x) { return num(x) && x + 0 >= 0 && x + 0 <= 1000000000 }
        function fail(message) {
            printf "%s:%d: %s\n", file, NR, message > "/dev/stderr"
            bad = 1
            exit 2
        }
        function normalized_comment(   text) {
            text = $0
            sub(/^#[[:space:]]*/, "", text)
            gsub(/[[:space:]]+/, " ", text)
            sub(/[[:space:]]+$/, "", text)
            return text
        }
        /^#/ {
            text = normalized_comment()
            if (text ~ /^schema=/) {
                if (text != "schema=whx-seedsearch-v1 primary_loss=0.10 secondary_loss=0.30 finalist_trials=paired tie_break=peel_then_lowest_seed")
                    fail("unsupported or mismatched seedsearch metadata")
                if (++schemas != 1) fail("duplicate seedsearch metadata")
            }
            else if (text ~ /^pairing=/) {
                count = split(text, fields, " ")
                trials = fields[2]
                sub(/^validation_trials=/, "", trials)
                if (count != 5 || fields[1] != "pairing=paired" ||
                    fields[2] !~ /^validation_trials=/ || !uint(trials) || trials + 0 < 1 ||
                    fields[3] != "trial_stride=2654435761" ||
                    fields[4] != "primary_base=0xBEE5+N*577" ||
                    fields[5] != "secondary_base=0xF00D+N*331")
                    fail("unsupported or malformed finalist pairing metadata")
                if (++pairings != 1) fail("duplicate finalist pairing metadata")
            }
            else if (text ~ /^N /) {
                if (text == "N best_pseed default_mean best_mean@0.10 best_mean@0.30") columns = 5
                else if (text == "N best_pseed best_dseed default_mean best_mean@0.10 best_mean@0.30") columns = 6
                else fail("unsupported seedsearch column header")
                if (++headers != 1) fail("duplicate seedsearch column header")
            }
            else if (text ~ /best_mean@/) {
                fail("mislabeled seedsearch loss metadata")
            }
            next
        }
        NF == 0 { next }
        {
            if (headers != 1) fail("data appears before the seedsearch column header")
            if (NF != columns) fail("data column count does not match header")
            if (!uint($1) || !uint($2)) fail("invalid N or peel seed")
            if ($1 + 0 < 2 || $1 + 0 > 64000 || $2 + 0 > 255) fail("N or peel seed out of range")
            if (columns == 5) {
                if (!score($3) || !score($4) || !score($5)) fail("invalid numeric score")
            }
            else {
                if (!sint($3) || $3 + 0 < -1 || $3 + 0 > 255 ||
                    !score($4) || !score($5) || !score($6)) fail("invalid dense seed or numeric score")
            }
            ++rows
        }
        END {
            if (bad) exit 2
            if (!bad && headers != 1) {
                printf "%s: missing seedsearch column header\n", file > "/dev/stderr"
                exit 2
            }
            if (!bad && rows == 0) {
                printf "%s: seedsearch input has no data rows\n", file > "/dev/stderr"
                exit 2
            }
            if (!bad && ((schemas == 1 && pairings != 1) || (schemas == 0 && pairings != 0))) {
                printf "%s: schema and finalist pairing metadata must appear together\n", file > "/dev/stderr"
                exit 2
            }
        }
    ' "$input"
done
