#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${OUT:-$ROOT/experiments/tiling/rowop_tiling}"
CXX="${CXX:-g++}"

"$CXX" -O3 -march=native -std=c++11 -Wall -Wextra -pthread \
    -I"$ROOT" \
    "$ROOT/experiments/tiling/rowop_tiling.cpp" \
    "$ROOT/gf256.cpp" \
    -o "$OUT"

echo "built $OUT"
