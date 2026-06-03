#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${OUT:-$ROOT/experiments/peeling/peel_sweep}"
CXX="${CXX:-g++}"

"$CXX" -O3 -march=native -std=c++11 -pthread -Wall -Wextra \
    -Wno-unused-function -Wno-implicit-fallthrough \
    -I"$ROOT" -I"$ROOT/include" \
    "$ROOT/experiments/peeling/peel_sweep.cpp" \
    "$ROOT/WirehairTools.cpp" \
    -o "$OUT"

echo "built $OUT"
