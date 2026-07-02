#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${OUT:-$ROOT/experiments/precode/precode_sim}"
CXX="${CXX:-g++}"

"$CXX" -O3 -march=native -std=c++11 -pthread -Wall -Wextra \
    -Wno-unused-function \
    -I"$ROOT" -I"$ROOT/include" \
    "$ROOT/experiments/precode/precode_sim.cpp" \
    "$ROOT/WirehairTools.cpp" \
    "$ROOT/codec/WirehairV2Precode.cpp" \
    "$ROOT/gf256.cpp" \
    -o "$OUT"

echo "built $OUT"
