cd /home/catid/wirehair
R=bench/results
while ps -p 182859 >/dev/null 2>&1; do sleep 20; done
echo "scan done $(date -u +%H:%M:%S)" > $R/task4.log
grep '^WEAK' $R/scan_17k_32k.log | sed -E 's/WEAK N=([0-9]+).*/\1/' | sort -n -u > $R/weakN_17k.txt
echo "weak N to seedsearch: $(wc -l < $R/weakN_17k.txt)" >> $R/task4.log
./bench/whx.knobfix seedsearch --threads 80 --nfile $R/weakN_17k.txt \
    --nseeds 32 --tsearch 150 --tverify 1800 > $R/seedsearch_17k.out 2>>$R/task4.log
echo "seedsearch done $(date -u +%H:%M:%S)" >> $R/task4.log
