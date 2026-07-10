#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
GEN_TABLES=${1:-"$ROOT/build/gen_tables"}

if [[ ! -x "$GEN_TABLES" ]]; then
    echo "missing gen_tables executable: $GEN_TABLES" >&2
    exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

expect_rejected()
{
    local name=$1
    shift

    if "$GEN_TABLES" --no-benchmarks --heavy-trials "$@" \
        > "$work/$name.out" 2> "$work/$name.err"
    then
        echo "gen_tables accepted invalid --heavy-trials value: $name" >&2
        exit 1
    fi
    if [[ -s "$work/$name.out" ]]; then
        echo "gen_tables started work before rejecting --heavy-trials: $name" >&2
        exit 1
    fi
    grep -Fq 'Usage:' "$work/$name.err"
}

expect_rejected missing
expect_rejected empty ''
expect_rejected negative -1
expect_rejected positive-sign +1
expect_rejected leading-space ' 1'
expect_rejected trailing-space '1 '
expect_rejected trailing-character 1x
expect_rejected hexadecimal 0x10
expect_rejected above-maximum 10000001
expect_rejected uint32-overflow 4294967296
expect_rejected uint64-overflow 18446744073709551616

# Parse the upper boundary without launching the deliberately large trial run.
"$GEN_TABLES" --heavy-trials 10000000 --help \
    > "$work/maximum.out" 2> "$work/maximum.err"
if [[ -s "$work/maximum.out" ]]; then
    echo "--help unexpectedly started the table generator" >&2
    exit 1
fi
grep -Fq 'maximum 10000000' "$work/maximum.err"

"$GEN_TABLES" --no-benchmarks --heavy-trials 4096 > "$work/first.out"
"$GEN_TABLES" --no-benchmarks --heavy-trials 4096 > "$work/second.out"
cmp "$work/first.out" "$work/second.out"

grep -Fqx 'Shipped Cauchy matrix reproduced from seed = 2318331135281' \
    "$work/first.out"
grep -Fqx '* Empirical perturbation singular rate: 14 / 4096 (expected ~1/256 = 16)' \
    "$work/first.out"
grep -Fqx '* Informational measurement only; not a reliability guarantee' \
    "$work/first.out"
if grep -Fq 'Tests failed' "$work/first.out"; then
    echo "gen_tables reported a failed regression gate" >&2
    exit 1
fi

echo "heavy-matrix generator regression passed (14/4096 empirical singular)"
