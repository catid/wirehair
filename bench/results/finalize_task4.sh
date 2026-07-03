cd /home/catid/wirehair
R=bench/results
while ! grep -q "seedsearch done" $R/task4.log 2>/dev/null; do sleep 30; done
echo "finalize start $(date -u +%H:%M:%S)" > $R/finalize4.log
# Extended table = original [2048,17500] + new [17500,32000] seedsearch outputs
DEFTHRESH=0 GOOD=0.045 GOOD30=0.15 OUT=WirehairPeelFixups.inc \
  bash bench/gen_fixups.sh $R/seedsearch.out $R/seedsearch_tail.out $R/seedsearch_17k.out >> $R/finalize4.log 2>&1
echo "table entries: $(grep -c '{ ' WirehairPeelFixups.inc) (range $(grep -oE '\{ *[0-9]+' WirehairPeelFixups.inc|grep -oE '[0-9]+'|sort -n|sed -n '1p;$p'|tr '\n' '-'))" >> $R/finalize4.log
OUT=bench/whx.fix TAG=fix bash bench/build.sh >> $R/finalize4.log 2>&1
# correctness gate
./bench/whx.fix fuzz --threads 70 --secs 35 --nmax 32000 --seed 0xF4 2>&1 | grep -E "fuzz done|fails=" | grep -v HIGH >> $R/finalize4.log
# spot-check a few previously-severe large-N defects now use a fixed seed (re-scan)
./bench/whx.fix scan --threads 70 --nlo 24080 --nhi 24420 --trials 200 --bb 64 --startmode 0 --loss 0.10 --thresh 0.10 2>/dev/null | grep -E "scan done|^WEAK" | head -8 >> $R/finalize4.log
echo "finalize done $(date -u +%H:%M:%S)" >> $R/finalize4.log
