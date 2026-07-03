cd /home/catid/wirehair
R=bench/results
echo "scan start $(date -u +%H:%M:%S)" > $R/ext.log
OUT=bench/whx.knobfix EXTRA=-DWH_SEED_KNOBS bash bench/build.sh >> $R/ext.log 2>&1
./bench/whx.knobfix scan --threads 80 --nlo 25787 --nhi 32000 --trials 150 --bb 64 --startmode 0 --loss 0.10 --thresh 0.10 > $R/scan_25k_32k.log 2>>$R/ext.log
echo "scan done $(date -u +%H:%M:%S), weak=$(grep -c '^WEAK' $R/scan_25k_32k.log)" >> $R/ext.log
grep '^WEAK' $R/scan_25k_32k.log | sed -E 's/WEAK N=([0-9]+).*/\1/' | sort -n -u > $R/weakN_25k.txt
./bench/whx.knobfix seedsearch --threads 80 --nfile $R/weakN_25k.txt --nseeds 40 --dseeds 32 --tsearch 150 --tverify 2000 --goodthr 0.04 > $R/joint_25k.out 2>>$R/ext.log
echo "joint done $(date -u +%H:%M:%S), $(grep -cE '^[0-9]' $R/joint_25k.out) N" >> $R/ext.log
GOOD=0.05 GOOD30=0.10 DEFTHRESH=0.05 bash bench/gen_fixups2.sh \
  $R/seedsearch.out $R/seedsearch_tail.out $R/seedsearch_17k.out $R/joint_hard.out $R/refine.out $R/joint_25k.out >> $R/ext.log 2>&1
OUT=bench/whx.fixed TAG=fixed bash bench/build.sh >> $R/ext.log 2>&1
./bench/whx.fixed fuzz --threads 70 --secs 30 --nmax 32000 --seed 0xC0DE 2>&1 | grep -E "overhead\(" | grep -v HIGH >> $R/ext.log
echo "ext done $(date -u +%H:%M:%S)" >> $R/ext.log
