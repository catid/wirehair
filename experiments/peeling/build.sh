#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${OUT:-$ROOT/experiments/peeling/peel_sweep}"
XOR_OUT="${XOR_OUT:-$ROOT/experiments/peeling/xor_bench}"
CXX="${CXX:-g++}"

"$CXX" -O3 -march=native -std=c++11 -pthread -Wall -Wextra \
    -Wno-unused-function -Wno-implicit-fallthrough \
    -I"$ROOT" -I"$ROOT/include" \
    "$ROOT/experiments/peeling/peel_sweep.cpp" \
    "$ROOT/WirehairTools.cpp" \
    -o "$OUT"

echo "built $OUT"

"$CXX" -O3 -march=native -std=c++11 -pthread -Wall -Wextra \
    "$ROOT/experiments/peeling/xor_bench.cpp" \
    -o "$XOR_OUT"

echo "built $XOR_OUT"
