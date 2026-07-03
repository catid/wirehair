cd /home/catid/wirehair
R=bench/results
echo "refine start $(date -u +%H:%M:%S)" > $R/refine.log
./bench/whx.knobfix seedsearch --threads 78 --nfile $R/refineN.txt \
    --nseeds 40 --dseeds 32 --tsearch 150 --tverify 2000 --goodthr 0.04 > $R/refine.out 2>>$R/refine.log
echo "refine search done $(date -u +%H:%M:%S), $(grep -cE '^[0-9]' $R/refine.out) N" >> $R/refine.log
GOOD=0.05 GOOD30=0.10 DEFTHRESH=0.05 bash bench/gen_fixups2.sh \
  $R/seedsearch.out $R/seedsearch_tail.out $R/seedsearch_17k.out $R/joint_hard.out $R/refine.out >> $R/refine.log 2>&1
OUT=bench/whx.fix TAG=fix bash bench/build.sh >> $R/refine.log 2>&1
./bench/whx.fix fuzz --threads 70 --secs 35 --nmax 32000 --seed 0xF6 2>&1 | grep -E "fuzz done|overhead\(" | grep -v HIGH >> $R/refine.log
for N in 24146 24166 6107 7533 11446 5550; do printf "  N=%s @0.30: " $N >> $R/refine.log
  ./bench/whx.fix ohead --threads 45 --trials 3000 --bb 64 --startmode 0 --loss 0.30 --nlo $N --nhi $N --nstep 1 2>&1|grep -vE '^#|^N '|awk '{print "mean="$2" max="$6}' >> $R/refine.log; done
echo "refine finalize done $(date -u +%H:%M:%S)" >> $R/refine.log
