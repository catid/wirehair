#!/bin/bash
# Verify a freshly-applied WirehairPeelFixups.inc:
#  1) rebuild,
#  2) re-scan the previously-weak N -> confirm overhead now normal (defects fixed by default),
#  3) correctness fuzz -> confirm decoding still exact (seed changes must not break recovery),
#  4) spot-check a band of non-fixed N -> confirm no collateral change.
set -e
cd "$(cd "$(dirname "$0")/.." && pwd)"
R=bench/results
echo "=== rebuild base with fixups ==="
OUT=bench/whx.fix TAG=fix bash bench/build.sh 2>&1 | grep -E "built|error:"

echo "=== (2) re-scan previously-weak N (should now report 0 WEAK) ==="
# weakN.txt has the original weak N; scan exactly those via tiny ranges is awkward, so
# scan the covered band again at the same threshold and compare counts.
./bench/whx.fix scan --threads 96 --nlo 2048 --nhi 24000 --trials 200 --bb 64 \
    --startmode 0 --loss 0.10 --thresh 0.10 > $R/scan_after.log 2>$R/scan_after.err
echo "weak BEFORE: $(grep -c '^WEAK' $R/scan_2k_24k.log)   weak AFTER: $(grep -c '^WEAK' $R/scan_after.log)"
echo "residual weak after fix:"; grep '^WEAK' $R/scan_after.log | head -20

echo "=== (3) correctness fuzz on fixed build (must be 0 fails) ==="
./bench/whx.fix fuzz --threads 96 --secs 30 --nmax 24000 --seed 0xF1XED 2>&1 | grep -vE "HIGH-OVERHEAD" | tail -3
