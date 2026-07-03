cd /home/catid/wirehair
R=bench/results
while ps -p 182859 >/dev/null 2>&1; do sleep 20; done   # exact PID, no pgrep self-match
echo "scan done $(date -u +%H:%M:%S)" > $R/task4.log
grep '^WEAK' $R/scan_17k_32k.log | sed -E 's/WEAK N=([0-9]+).*/\1/' | sort -n -u > $R/weakN_17k.txt
echo "weak N: $(wc -l < $R/weakN_17k.txt)" >> $R/task4.log
# seedsearch (uncapped/normal distribution); fast params for large N
./bench/whx.knobfix seedsearch --threads 80 --nfile $R/weakN_17k.txt \
    --nseeds 40 --tsearch 200 --tverify 2500 > $R/seedsearch_17k.out 2>>$R/task4.log
echo "seedsearch done $(date -u +%H:%M:%S)" >> $R/task4.log
