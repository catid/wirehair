#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${OUT:-$ROOT/experiments/tiling/rowop_tiling}"
CXX="${CXX:-g++}"

"$CXX" -O3 -march=native -std=c++11 -Wall -Wextra -pthread \
    "$ROOT/experiments/tiling/rowop_tiling.cpp" \
    -o "$OUT"

echo "built $OUT"
