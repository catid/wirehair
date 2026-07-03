#!/bin/bash
set -u
set -o pipefail

cd /home/catid/wirehair
FAIL_FILE=/tmp/queue_fail
rm -f "$FAIL_FILE" /tmp/e5small_done /tmp/e3abig_done /tmp/e1bbig_done /tmp/e5big_done /tmp/queue_done

mark_done() {
  local marker="$1"
  if [ -s "$FAIL_FILE" ]; then
    echo "queue failures before $marker:" >&2
    cat "$FAIL_FILE" >&2
    exit 1
  fi
  echo OK > "$marker"
}

./experiments/precode/precode_sim --K 1000 --oh 0,1,2 --trials 20000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s50_d12,ldpcdense_s50_d12_s2,ldpcdense_n12_s50_d12,ldpcdense_n12_s50_d12_s2,ldpcdense_n14_s50_d12,ldpcdense_s50_d15 \
  > experiments/precode/results/e5_K1000.csv 2>/tmp/q_e5_K1000.csv.err || echo "FAIL e5_K1000.csv" >> /tmp/queue_fail
./experiments/precode/precode_sim --K 3200 --oh 0,1,2 --trials 20000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s62_d12,ldpcdense_s62_d12_s2,ldpcdense_n12_s62_d12,ldpcdense_n12_s62_d12_s2,ldpcdense_n14_s62_d12,ldpcdense_s62_d19 \
  > experiments/precode/results/e5_K3200.csv 2>/tmp/q_e5_K3200.csv.err || echo "FAIL e5_K3200.csv" >> /tmp/queue_fail
./experiments/precode/precode_sim --K 10000 --oh 0,1,2 --trials 20000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s86_d12,ldpcdense_s86_d12_s2,ldpcdense_n12_s86_d12,ldpcdense_n12_s86_d12_s2,ldpcdense_n14_s86_d12,ldpcdense_s86_d26 \
  > experiments/precode/results/e5_K10000.csv 2>/tmp/q_e5_K10000.csv.err || echo "FAIL e5_K10000.csv" >> /tmp/queue_fail
mark_done /tmp/e5small_done
./experiments/precode/precode_sim --K 10000 --oh 0,1,2,5 --trials 8000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s52_d12,ldpcdense_s52_d17,ldpcdense_s52_d26,ldpcdense_s60_d12,ldpcdense_s60_d17,ldpcdense_s60_d26,ldpcdense_s69_d12,ldpcdense_s69_d17,ldpcdense_s69_d26,ldpcdense_s77_d12,ldpcdense_s77_d17,ldpcdense_s77_d26,ldpcdense_s86_d12,ldpcdense_s86_d17,ldpcdense_s86_d26,ldpcdense_s69_d17_h12,ldpcdense_s86_d12_h12 \
  > experiments/precode/results/e3a_K10000.csv 2>/tmp/q_e3a_K10000.csv.err || echo "FAIL e3a_K10000.csv" >> /tmp/queue_fail
./experiments/precode/precode_sim --K 16000 --oh 0,1,2,5 --trials 8000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s68_d12,ldpcdense_s68_d23,ldpcdense_s68_d34,ldpcdense_s80_d12,ldpcdense_s80_d23,ldpcdense_s80_d34,ldpcdense_s91_d12,ldpcdense_s91_d23,ldpcdense_s91_d34,ldpcdense_s103_d12,ldpcdense_s103_d23,ldpcdense_s103_d34,ldpcdense_s114_d12,ldpcdense_s114_d23,ldpcdense_s114_d34,ldpcdense_s91_d23_h12,ldpcdense_s114_d12_h12 \
  > experiments/precode/results/e3a_K16000.csv 2>/tmp/q_e3a_K16000.csv.err || echo "FAIL e3a_K16000.csv" >> /tmp/queue_fail
./experiments/precode/precode_sim --K 20000 --oh 0,1,2,5 --trials 8000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s80_d12,ldpcdense_s80_d27,ldpcdense_s80_d40,ldpcdense_s94_d12,ldpcdense_s94_d27,ldpcdense_s94_d40,ldpcdense_s107_d12,ldpcdense_s107_d27,ldpcdense_s107_d40,ldpcdense_s121_d12,ldpcdense_s121_d27,ldpcdense_s121_d40,ldpcdense_s134_d12,ldpcdense_s134_d27,ldpcdense_s134_d40,ldpcdense_s107_d27_h12,ldpcdense_s134_d12_h12 \
  > experiments/precode/results/e3a_K20000.csv 2>/tmp/q_e3a_K20000.csv.err || echo "FAIL e3a_K20000.csv" >> /tmp/queue_fail
./experiments/precode/precode_sim --K 32000 --oh 0,1,2,5 --trials 5000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s114_d12,ldpcdense_s114_d38,ldpcdense_s114_d57,ldpcdense_s133_d12,ldpcdense_s133_d38,ldpcdense_s133_d57,ldpcdense_s152_d12,ldpcdense_s152_d38,ldpcdense_s152_d57,ldpcdense_s171_d12,ldpcdense_s171_d38,ldpcdense_s171_d57,ldpcdense_s190_d12,ldpcdense_s190_d38,ldpcdense_s190_d57,ldpcdense_s152_d38_h12,ldpcdense_s190_d12_h12 \
  > experiments/precode/results/e3a_K32000.csv 2>/tmp/q_e3a_K32000.csv.err || echo "FAIL e3a_K32000.csv" >> /tmp/queue_fail
./experiments/precode/precode_sim --K 48000 --oh 0,1,2,5 --trials 3000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s222_d12,ldpcdense_s222_d74,ldpcdense_s222_d111,ldpcdense_s259_d12,ldpcdense_s259_d74,ldpcdense_s259_d111,ldpcdense_s296_d12,ldpcdense_s296_d74,ldpcdense_s296_d111,ldpcdense_s333_d12,ldpcdense_s333_d74,ldpcdense_s333_d111,ldpcdense_s370_d12,ldpcdense_s370_d74,ldpcdense_s370_d111,ldpcdense_s296_d74_h12,ldpcdense_s370_d12_h12 \
  > experiments/precode/results/e3a_K48000.csv 2>/tmp/q_e3a_K48000.csv.err || echo "FAIL e3a_K48000.csv" >> /tmp/queue_fail
./experiments/precode/precode_sim --K 64000 --oh 0,1,2,5 --trials 3000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s208_d12,ldpcdense_s208_d69,ldpcdense_s208_d104,ldpcdense_s242_d12,ldpcdense_s242_d69,ldpcdense_s242_d104,ldpcdense_s277_d12,ldpcdense_s277_d69,ldpcdense_s277_d104,ldpcdense_s311_d12,ldpcdense_s311_d69,ldpcdense_s311_d104,ldpcdense_s346_d12,ldpcdense_s346_d69,ldpcdense_s346_d104,ldpcdense_s277_d69_h12,ldpcdense_s346_d12_h12 \
  > experiments/precode/results/e3a_K64000.csv 2>/tmp/q_e3a_K64000.csv.err || echo "FAIL e3a_K64000.csv" >> /tmp/queue_fail
mark_done /tmp/e3abig_done
./experiments/precode/precode_sim --K 32000 --oh 0,1,2 --trials 3000 --threads 116 \
  --paired --seed 24414209 --ge-replay \
  --schemes dense,dense_h12,ldpc,ldpc_h12,ldpcdense_s152_d38,ldpcdense_s152_d38_h8,ldpcdense_s152_d38_h12,ldpcdense_s190_d12_h12 \
  > experiments/precode/results/e1_hgrid_K32000.csv 2>/tmp/q_e1_hgrid_K32000.csv.err || echo "FAIL e1_hgrid_K32000.csv" >> /tmp/queue_fail
./experiments/precode/precode_sim --K 64000 --oh 0,1,2 --trials 1500 --threads 116 \
  --paired --seed 24414209 --ge-replay \
  --schemes dense,dense_h12,ldpc,ldpc_h12,ldpcdense_s277_d69,ldpcdense_s277_d69_h8,ldpcdense_s277_d69_h12,ldpcdense_s346_d12_h12 \
  > experiments/precode/results/e1_hgrid_K64000.csv 2>/tmp/q_e1_hgrid_K64000.csv.err || echo "FAIL e1_hgrid_K64000.csv" >> /tmp/queue_fail
mark_done /tmp/e1bbig_done
./experiments/precode/precode_sim --K 32000 --oh 0,1,2 --trials 5000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpcdense_s190_d12,ldpcdense_s190_d12_s2,ldpcdense_n12_s190_d12,ldpcdense_n12_s190_d12_s2,ldpcdense_s190_d57 \
  > experiments/precode/results/e5_K32000.csv 2>/tmp/q_e5_K32000.csv.err || echo "FAIL e5_K32000.csv" >> /tmp/queue_fail
./experiments/precode/precode_sim --K 64000 --oh 0,1,2 --trials 3000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpcdense_s346_d12,ldpcdense_s346_d12_s2,ldpcdense_n12_s346_d12,ldpcdense_n12_s346_d12_s2,ldpcdense_s346_d104 \
  > experiments/precode/results/e5_K64000.csv 2>/tmp/q_e5_K64000.csv.err || echo "FAIL e5_K64000.csv" >> /tmp/queue_fail
mark_done /tmp/e5big_done
mark_done /tmp/queue_done
