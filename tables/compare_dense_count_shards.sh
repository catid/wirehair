#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CXX=${CXX:-g++}

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

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

"$CXX" -std=c++11 -O2 \
    -I"$ROOT_DIR" \
    -I"$ROOT_DIR/include" \
    -o "$tmpdir/compare_dense_count" \
    -x c++ - \
    "$ROOT_DIR/WirehairTools.cpp" \
    "$ROOT_DIR/gf256.cpp" <<'CPP'
#include "WirehairTools.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

int main()
{
    std::string line;
    unsigned total = 0;
    unsigned changed = 0;
    unsigned lower = 0;
    unsigned higher = 0;
    int min_delta = 0;
    int max_delta = 0;
    long long delta_sum = 0;

    std::cout << "Source\tN\tCurrentDense\tGeneratedDense\tDelta\n";
    while (std::getline(std::cin, line))
    {
        const std::string::size_type tab1 = line.find('\t');
        const std::string::size_type tab2 =
            tab1 == std::string::npos ? std::string::npos :
            line.find('\t', tab1 + 1);
        if (tab1 == std::string::npos || tab2 == std::string::npos) {
            continue;
        }
        const std::string source = line.substr(0, tab1);
        const unsigned n = (unsigned)std::strtoul(
            line.c_str() + tab1 + 1, 0, 10);
        const unsigned generated = (unsigned)std::strtoul(
            line.c_str() + tab2 + 1, 0, 10);
        const unsigned current = wirehair::GetDenseCount(n);
        const int delta = (int)generated - (int)current;
        ++total;
        if (delta == 0) {
            continue;
        }
        if (changed == 0 || delta < min_delta) {
            min_delta = delta;
        }
        if (changed == 0 || delta > max_delta) {
            max_delta = delta;
        }
        ++changed;
        if (delta < 0) {
            ++lower;
        }
        else {
            ++higher;
        }
        delta_sum += delta;
        std::cout << source << '\t' << n << '\t' << current << '\t'
                  << generated << '\t' << delta << '\n';
    }

    if (total == 0) {
        std::cerr << "no dense-count rows found\n";
        return 1;
    }

    const double mean_delta =
        changed == 0 ? 0.0 : (double)delta_sum / (double)changed;
    std::cerr << "total=" << total
              << " changed=" << changed
              << " lower=" << lower
              << " higher=" << higher
              << " min_delta=" << min_delta
              << " max_delta=" << max_delta
              << " mean_changed_delta=" << std::fixed << std::setprecision(2)
              << mean_delta << '\n';
    return 0;
}
CPP

for input in "${inputs[@]}"; do
    awk -v source="$input" '
        NR == 1 { next }
        NF >= 2 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ {
            printf "%s\t%s\t%s\n", source, $1, $2
        }
    ' "$input"
done | "$tmpdir/compare_dense_count"
