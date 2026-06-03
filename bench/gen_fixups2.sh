#!/bin/bash
# Generate WirehairPeelFixups.inc + WirehairDenseFixups.inc from seedsearch outputs.
# Formats: 5-col "N pseed def b10 b30" (peel-only), 6-col "N pseed dseed def b10 b30" (joint).
# When a N appears in multiple outputs, keep the BEST (peel,dense) combo (min of max(b10,b30)).
# Keep only genuine verified fixes (default weak; best good at loss 0.10 AND 0.30). Sorted by N.
set -e
GOOD=${GOOD:-0.05}
GOOD30=${GOOD30:-0.10}
DEFTHRESH=${DEFTHRESH:-0.05}
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
TMP="$(mktemp)"

cat "$@" | awk -v dt="$DEFTHRESH" -v g="$GOOD" -v g3="$GOOD30" '
  function mx(a,b){ return a>b?a:b }
  /^#/ { next }
  NF==5 { N=$1; p=$2; d=-1; def=$3; b10=$4; b30=$5 }
  NF==6 { N=$1; p=$2; d=$3; def=$4; b10=$5; b30=$6 }
  NF==5 || NF==6 {
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
} > WirehairPeelFixups.inc

{
  echo "// Dense-seed corrections (joint peel+dense search, Task5). Paired with the peel"
  echo "// correction for the same N. Sorted ascending by N (binary search). Format: { N, seed }"
  awk '$3>=0 {printf "    { %5d, %3d },\n", $1, $3}' "$TMP"
} > WirehairDenseFixups.inc

grep -q '^    {' WirehairPeelFixups.inc  || echo "    { 1, 0 }," >> WirehairPeelFixups.inc
grep -q '^    {' WirehairDenseFixups.inc || echo "    { 1, 0 }," >> WirehairDenseFixups.inc
rm -f "$TMP"
echo "peel corrections: $(grep -c '^    {' WirehairPeelFixups.inc) | dense corrections: $(grep -c '^    {' WirehairDenseFixups.inc)"
