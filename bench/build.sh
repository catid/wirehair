#!/bin/bash
# Build the whx harness by compiling wirehair sources directly.
# Env:
#   EXTRA  : extra compiler flags / defines (e.g. "-DWH_GFNI=1")
#   LINK_EXTRA : extra linker flags
#   OUT    : output binary path (default bench/whx)
#   TAG    : object subdir tag for parallel A/B builds (default "base")
#   CXX    : C++ compiler (default g++)
#   LTO    : 0/off, 1/on/full, auto, or thin
#   PGO    : off, gen/generate, or use
#   PGO_DIR: profile directory (default bench/pgo/$TAG)
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT="${OUT:-bench/whx}"
TAG="${TAG:-base}"
OBJ="bench/obj/$TAG"
CXX="${CXX:-g++}"
LTO="${LTO:-0}"
PGO="${PGO:-off}"
PGO_DIR="${PGO_DIR:-bench/pgo/$TAG}"
mkdir -p "$OBJ"
FLAGS="-O3 -march=native -std=c++11 -pthread -Wall -Wno-unused-function -I$ROOT -I$ROOT/include -I$ROOT/test"
LINK_FLAGS=""

case "$LTO" in
    0|off|OFF|false|FALSE|"") ;;
    1|on|ON|true|TRUE|full|FULL)
        FLAGS="$FLAGS -flto"
        LINK_FLAGS="$LINK_FLAGS -flto"
        ;;
    auto|AUTO)
        FLAGS="$FLAGS -flto=auto"
        LINK_FLAGS="$LINK_FLAGS -flto=auto"
        ;;
    thin|THIN)
        if ! "$CXX" -dM -E -x c++ /dev/null 2>/dev/null | grep -q '__clang__'; then
            echo "LTO=thin requires Clang; CXX='$CXX' does not support -flto=thin" >&2
            exit 2
        fi
        FLAGS="$FLAGS -flto=thin"
        LINK_FLAGS="$LINK_FLAGS -flto=thin"
        ;;
    *)
        echo "Unknown LTO='$LTO' (expected off/on/auto/thin)" >&2
        exit 2
        ;;
esac

case "$PGO" in
    off|OFF|0|"") ;;
    gen|GEN|generate|GENERATE)
        mkdir -p "$PGO_DIR"
        FLAGS="$FLAGS -fprofile-generate=$PGO_DIR"
        LINK_FLAGS="$LINK_FLAGS -fprofile-generate=$PGO_DIR"
        ;;
    use|USE)
        FLAGS="$FLAGS -fprofile-use=$PGO_DIR -fprofile-correction -Wno-missing-profile"
        LINK_FLAGS="$LINK_FLAGS -fprofile-use=$PGO_DIR -fprofile-correction -Wno-missing-profile"
        ;;
    *)
        echo "Unknown PGO='$PGO' (expected off/gen/use)" >&2
        exit 2
        ;;
esac

FLAGS="$FLAGS $EXTRA"
LINK_FLAGS="$LINK_FLAGS $LINK_EXTRA"

pids=""
for src in gf256 WirehairCodec WirehairTools wirehair; do
    "$CXX" $FLAGS -c "$src.cpp" -o "$OBJ/$src.o" & pids="$pids $!"
done
"$CXX" $FLAGS -c bench/whx.cpp -o "$OBJ/whx.o" & pids="$pids $!"
fail=0
for p in $pids; do wait "$p" || fail=1; done
if [ "$fail" != "0" ]; then echo "COMPILE FAILED"; exit 1; fi
"$CXX" "$OBJ"/gf256.o "$OBJ"/WirehairCodec.o "$OBJ"/WirehairTools.o "$OBJ"/wirehair.o "$OBJ"/whx.o -pthread $LINK_FLAGS -o "$OUT"
echo "built $OUT  (CXX='$CXX' EXTRA='$EXTRA' LTO='$LTO' PGO='$PGO' PGO_DIR='$PGO_DIR')"
