#!/bin/bash
# Chains after the scan: waits for it to finish (one campaign at a time), extracts weak N,
# runs a small validation, then the full seedsearch. Writes to bench/results/.
set -euo pipefail
cd "$(cd "$(dirname "$0")/.." && pwd)"
R=bench/results
# Wait for the scan to appear before waiting for it to finish: if it never
# started (or already crashed) this used to fall through instantly and
# process a partial scan log as if complete.
for _ in $(seq 60); do
    pgrep -f "whx.base scan" >/dev/null && break
    sleep 1
done
pgrep -f "whx.base scan" >/dev/null || { echo "scan not running; aborting" >&2; exit 1; }
while pgrep -f "whx.base scan" >/dev/null; do sleep 15; done
echo "scan finished $(date -u +%H:%M:%S)" > $R/drive.log
grep '^WEAK' $R/scan_2k_24k.log | sed -E 's/WEAK N=([0-9]+).*/\1/' | sort -n -u > $R/weakN.txt || true
echo "weak N count: $(wc -l < $R/weakN.txt)" >> $R/drive.log
if [ ! -s "$R/weakN.txt" ]; then
    echo "NO_WEAK_DONE $(date -u +%H:%M:%S)" >> $R/drive.log
    exit 0
fi
# validation: weak {2962,3061,3430} + good {3060,3062,5000} -- sanity-check parallel verify
./bench/whx.knobfix seedsearch --threads 96 --nlist "2962,3061,3430,3060,3062,5000" \
    --nseeds 96 --tsearch 400 --tverify 5000 > $R/seedsearch_val.out 2>>$R/drive.log
echo "VAL_DONE $(date -u +%H:%M:%S)" >> $R/drive.log
# full campaign over all weak N (lighter params: almost any non-default seed works, so a
# modest search + robust parallel verify suffices)
./bench/whx.knobfix seedsearch --threads 96 --nfile $R/weakN.txt \
    --nseeds 96 --tsearch 350 --tverify 5000 > $R/seedsearch.out 2>>$R/drive.log
echo "FULL_DONE $(date -u +%H:%M:%S)" >> $R/drive.log
