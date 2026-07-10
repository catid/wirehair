#!/bin/bash
# Generate WirehairPeelFixups.inc + WirehairDenseFixups.inc from seedsearch outputs.
# Formats: 5-col "N pseed def b10 b30" (peel-only), 6-col "N pseed dseed def b10 b30" (joint).
# When a N appears in multiple outputs, keep the BEST (peel,dense) combo (min of max(b10,b30)).
# Keep only genuine verified fixes (default weak; best good at loss 0.10 AND 0.30). Sorted by N.
set -euo pipefail
GOOD=${GOOD:-0.05}
GOOD30=${GOOD30:-0.10}
DEFTHRESH=${DEFTHRESH:-0.05}
PEEL_OUT=${PEEL_OUT:-WirehairPeelFixups.inc}
DENSE_OUT=${DENSE_OUT:-WirehairDenseFixups.inc}
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
bash bench/validate_seedsearch.sh "$@"
TMP="$(mktemp)"
PEEL_TMP="$(mktemp "${PEEL_OUT}.tmp.XXXXXX")"
DENSE_TMP="$(mktemp "${DENSE_OUT}.tmp.XXXXXX")"
trap 'rm -f "$TMP" "$PEEL_TMP" "$DENSE_TMP"' EXIT

cat "$@" | awk -v dt="$DEFTHRESH" -v g="$GOOD" -v g3="$GOOD30" '
  function uint(x){ return x ~ /^[0-9]+$/ }
  function sint(x){ return x ~ /^-?[0-9]+$/ }
  function num(x){ return x ~ /^([0-9]+([.][0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?$/ }
  function die(){ printf "malformed seedsearch line: %s\n", $0 > "/dev/stderr"; exit 2 }
  function mx(a,b){ return a>b?a:b }
  /^#/ || NF==0 { next }
  NF!=5 && NF!=6 { die() }
  NF==5 {
    if (!uint($1) || !uint($2) || !num($3) || !num($4) || !num($5)) die();
    N=$1; p=$2; d=-1; def=$3; b10=$4; b30=$5
  }
  NF==6 {
    if (!uint($1) || !uint($2) || !sint($3) || !num($4) || !num($5) || !num($6)) die();
    N=$1; p=$2; d=$3; def=$4; b10=$5; b30=$6
  }
  NF==5 || NF==6 {
    if (N+0 < 2 || N+0 > 64000 || p+0 > 255 || d+0 < -1 || d+0 > 255) die();
    if (def+0>dt && b10+0<g && b30+0<g3 && N+0>=2 && N+0<=64000 && p+0>=0 && p+0<=255) {
      s = mx(b10+0, b30+0)
      if (!(N in bs) || s < bs[N]) { bs[N]=s; bp[N]=p+0; bd[N]=d+0 }
    }
  }
  END { for (n in bs) print n+0, bp[n], bd[n] }
' | sort -n -k1,1 > "$TMP"

{
  echo "// Peel-seed corrections (bench/whx seedsearch + gen_fixups2.sh; best combo per N)."
  echo "// Sorted ascending by N (binary search). Format: { N, seed }"
  awk '{printf "    { %5d, %3d },\n", $1, $2}' "$TMP"
} > "$PEEL_TMP"

{
  echo "// Dense-seed corrections (joint peel+dense search, Task5). Paired with the peel"
  echo "// correction for the same N. Sorted ascending by N (binary search). Format: { N, seed }"
  awk '$3>=0 {printf "    { %5d, %3d },\n", $1, $3}' "$TMP"
} > "$DENSE_TMP"

grep -q '^    {' "$PEEL_TMP" || echo "    { 1, 0 }," >> "$PEEL_TMP"
grep -q '^    {' "$DENSE_TMP" || echo "    { 1, 0 }," >> "$DENSE_TMP"
mv -f "$PEEL_TMP" "$PEEL_OUT"
mv -f "$DENSE_TMP" "$DENSE_OUT"
rm -f "$TMP"
trap - EXIT
echo "peel corrections: $(grep -c '^    {' "$PEEL_OUT") | dense corrections: $(grep -c '^    {' "$DENSE_OUT")"
