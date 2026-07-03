cd /home/catid/wirehair
while ps -p 121876 >/dev/null 2>&1; do sleep 15; done   # exact PID, no pgrep self-match
DEFTHRESH=0.05 GOOD=0.045 GOOD30=0.06 OUT=WirehairPeelFixups.inc \
  bash bench/gen_fixups.sh bench/results/seedsearch.out bench/results/seedsearch_tail.out > bench/results/finalize.log 2>&1
OUT=bench/whx.fix TAG=fix bash bench/build.sh >> bench/results/finalize.log 2>&1
echo "FINALIZER_DONE entries=$(grep -c '{ ' WirehairPeelFixups.inc)" >> bench/results/finalize.log
