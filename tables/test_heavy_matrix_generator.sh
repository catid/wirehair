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
