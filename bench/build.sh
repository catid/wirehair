#!/bin/bash
# Build the whx harness by compiling wirehair sources directly.
# Env:
#   EXTRA  : extra compiler flags / defines (e.g. "-DWH_GFNI=1")
#   OUT    : output binary path (default bench/whx)
#   TAG    : object subdir tag for parallel A/B builds (default "base")
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT="${OUT:-bench/whx}"
TAG="${TAG:-base}"
OBJ="bench/obj/$TAG"
mkdir -p "$OBJ"
FLAGS="-O3 -march=native -std=c++11 -pthread -Wall -Wno-unused-function -I$ROOT -I$ROOT/include -I$ROOT/test $EXTRA"

pids=""
for src in gf256 WirehairCodec WirehairTools wirehair; do
    g++ $FLAGS -c "$src.cpp" -o "$OBJ/$src.o" & pids="$pids $!"
done
g++ $FLAGS -c bench/whx.cpp -o "$OBJ/whx.o" & pids="$pids $!"
fail=0
for p in $pids; do wait "$p" || fail=1; done
if [ "$fail" != "0" ]; then echo "COMPILE FAILED"; exit 1; fi
g++ "$OBJ"/gf256.o "$OBJ"/WirehairCodec.o "$OBJ"/WirehairTools.o "$OBJ"/wirehair.o "$OBJ"/whx.o -pthread -o "$OUT"
echo "built $OUT  (EXTRA='$EXTRA')"
