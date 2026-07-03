cd /home/catid/wirehair
R=bench/results
while ps -p 629687 >/dev/null 2>&1; do sleep 20; done
echo "joint done $(date -u +%H:%M:%S)" > $R/finalize5.log
GOOD=0.05 GOOD30=0.15 DEFTHRESH=0.05 bash bench/gen_fixups2.sh \
  $R/seedsearch.out $R/seedsearch_tail.out $R/seedsearch_17k.out $R/joint_hard.out >> $R/finalize5.log 2>&1
OUT=bench/whx.fix TAG=fix bash bench/build.sh >> $R/finalize5.log 2>&1
./bench/whx.fix fuzz --threads 70 --secs 35 --nmax 32000 --seed 0xF5 2>&1 | grep -E "fuzz done|overhead\(" | grep -v HIGH >> $R/finalize5.log
for N in 7533 11446 24146; do
  printf "  N=%s @0.30: " $N >> $R/finalize5.log
  ./bench/whx.fix ohead --threads 45 --trials 3000 --bb 64 --startmode 0 --loss 0.30 --nlo $N --nhi $N --nstep 1 2>&1 | grep -vE '^#|^N ' | awk '{print "mean="$2" max="$6}' >> $R/finalize5.log
done
echo "finalize5 done $(date -u +%H:%M:%S)" >> $R/finalize5.log
