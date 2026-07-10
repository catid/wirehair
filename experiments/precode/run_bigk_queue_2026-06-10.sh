#!/bin/bash
set -u
set -o pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT" || exit 1
FAIL_FILE=/tmp/queue_fail
RESULT_DIR=experiments/precode/results
SIM=./experiments/precode/precode_sim
VALIDATOR=./experiments/precode/validate_results.py

build_sim() {
  if ! OUT="$SIM" bash experiments/precode/build.sh; then
    echo "FAIL build precode_sim" >> "$FAIL_FILE"
    cat "$FAIL_FILE" >&2
    exit 1
  fi
  if [ ! -x "$SIM" ] || [ "$SIM" -ot experiments/precode/precode_sim.cpp ]; then
    echo "FAIL precode_sim missing or stale after build" >> "$FAIL_FILE"
    cat "$FAIL_FILE" >&2
    exit 1
  fi
}

prepare_queue() {
  rm -f "$FAIL_FILE" /tmp/e5small_done /tmp/e3abig_done /tmp/e1bbig_done /tmp/e5big_done /tmp/queue_done
  mkdir -p "$RESULT_DIR"
  build_sim
}

run_csv() {
  local out="$1"
  shift
  local err="${out}.err"
  local tmp_csv
  local tmp_err
  tmp_csv="$(mktemp "${out}.tmp.XXXXXX")" || exit 1
  tmp_err="$(mktemp "${err}.tmp.XXXXXX")" || {
    rm -f "$tmp_csv"
    exit 1
  }

  "$@" > "$tmp_csv" 2> "$tmp_err"
  local rc=$?
  if [ "$rc" -eq 0 ]; then
    local expect_k="" expect_oh="" expect_schemes="" expect_trials=""
    local previous=""
    for argument in "$@"; do
      case "$previous" in
        --K) expect_k="$argument" ;;
        --oh) expect_oh="$argument" ;;
        --schemes) expect_schemes="$argument" ;;
        --trials) expect_trials="$argument" ;;
      esac
      previous="$argument"
    done
    local validation=(python3 "$VALIDATOR" "$tmp_csv")
    [ -n "$expect_k" ] && validation+=(--expect-k "$expect_k")
    [ -n "$expect_oh" ] && validation+=(--expect-oh "$expect_oh")
    [ -n "$expect_schemes" ] && validation+=(--expect-schemes "$expect_schemes")
    [ -n "$expect_trials" ] && validation+=(--expect-trials "$expect_trials")
    # Keep command-specific grid expectations on the candidate alone.  Then
    # reject cell collisions with final shards sharing its <campaign>_K prefix.
    local out_dir="."
    if [[ "$out" == */* ]]; then
      out_dir="${out%/*}"
    fi
    local out_name="${out##*/}"
    local campaign="${out_name%%_K*}"
    local campaign_validation=(python3 "$VALIDATOR")
    local shard
    while IFS= read -r -d '' shard; do
      if [ "$shard" != "$out" ]; then
        campaign_validation+=("$shard")
      fi
    done < <(find "$out_dir" -maxdepth 1 -type f \
      -name "${campaign}_K*.csv" -print0)
    campaign_validation+=("$tmp_csv")
    if ! "${validation[@]}" > /dev/null 2>> "$tmp_err" ||
        ! "${campaign_validation[@]}" > /dev/null 2>> "$tmp_err"; then
      mv -f "$tmp_csv" "${out}.invalid"
      mv -f "$tmp_err" "$err"
      echo "FAIL $out validation" >> "$FAIL_FILE"
      return 1
    fi
    if [ -s "$tmp_err" ]; then
      mv -f "$tmp_err" "$err"
      rm -f "$tmp_csv"
      echo "FAIL $out stderr non-empty" >> "$FAIL_FILE"
      return 1
    fi
    mv -f "$tmp_csv" "$out"
    rm -f "$tmp_err" "$err"
    return 0
  fi

  if [ -s "$tmp_err" ]; then
    mv -f "$tmp_err" "$err"
  else
    printf "command failed with exit %d\n" "$rc" > "$err"
    rm -f "$tmp_err"
  fi
  rm -f "$tmp_csv"
  echo "FAIL $out rc=$rc" >> "$FAIL_FILE"
  return 1
}

mark_done() {
  local marker="$1"
  if [ -s "$FAIL_FILE" ]; then
    echo "queue failures before $marker:" >&2
    cat "$FAIL_FILE" >&2
    exit 1
  fi
  echo OK > "$marker"
}

queue_self_test() {
  local test_dir
  test_dir="$(mktemp -d)" || exit 1
  # Capture the local path now; it is out of scope when the EXIT trap runs.
  # shellcheck disable=SC2064
  trap "rm -rf -- '$test_dir'" EXIT
  FAIL_FILE="$test_dir/failures"

  local valid="$test_dir/valid.csv"
  if ! run_csv "$valid" "$SIM" --K 32 --oh 0,1 --trials 2 --threads 1 \
      --paired --seed 123 --schemes dense; then
    echo "queue self-test: valid result was not promoted" >&2
    exit 1
  fi
  if [ ! -s "$valid" ] || [ -e "${valid}.err" ] ||
      find "$test_dir" -type f -name '*.tmp.*' -print -quit | grep -q .; then
    echo "queue self-test: valid promotion left missing or stale artifacts" >&2
    exit 1
  fi

  local existing="$test_dir/collision_Kexisting.csv"
  local collision="$test_dir/collision_Kcandidate.csv"
  cp -f "$valid" "$existing"
  if run_csv "$collision" sed -n p "$existing"; then
    echo "queue self-test: duplicate campaign cell was promoted" >&2
    exit 1
  fi
  if [ -e "$collision" ] || [ ! -s "${collision}.invalid" ] ||
      [ ! -s "${collision}.err" ] ||
      ! grep -Fq "first seen at $existing:2" "${collision}.err" ||
      ! grep -Eq 'collision_Kcandidate\.csv\.tmp\.[^:]+:2: duplicate experiment row' \
        "${collision}.err"; then
    echo "queue self-test: duplicate campaign cell was not quarantined" >&2
    exit 1
  fi

  local unrelated="$test_dir/unrelated_Kexisting.csv"
  local isolated="$test_dir/isolated_Kcandidate.csv"
  cp -f "$valid" "$unrelated"
  if ! run_csv "$isolated" sed -n p "$unrelated" || [ ! -s "$isolated" ] ||
      [ -e "${isolated}.err" ]; then
    echo "queue self-test: unrelated campaign blocked valid promotion" >&2
    exit 1
  fi

  local invalid="$test_dir/invalid.csv"
  if run_csv "$invalid" bash -c 'printf "K,scheme\n32,dense\nBROKEN\n"' \
      queue-fixture --K 32 --oh 0 --trials 2 --schemes dense; then
    echo "queue self-test: malformed result was promoted" >&2
    exit 1
  fi
  if [ -e "$invalid" ] || [ ! -s "${invalid}.invalid" ] ||
      [ ! -s "${invalid}.err" ] || [ ! -s "$FAIL_FILE" ]; then
    echo "queue self-test: malformed result was not quarantined" >&2
    exit 1
  fi

  local failed="$test_dir/failed.csv"
  if run_csv "$failed" bash -c 'exit 37' queue-fixture; then
    echo "queue self-test: failing command was promoted" >&2
    exit 1
  fi
  if [ -e "$failed" ] || [ -e "${failed}.invalid" ] ||
      ! grep -q 'rc=37' "$FAIL_FILE" ||
      ! grep -q 'exit 37' "${failed}.err"; then
    echo "queue self-test: command failure status was not preserved" >&2
    exit 1
  fi
  echo "queue result promotion self-test passed"
}

if [ "${1:-}" = "--self-test" ]; then
  build_sim
  queue_self_test
  exit 0
fi
prepare_queue

run_csv "$RESULT_DIR/e5_K1000.csv" "$SIM" --K 1000 --oh 0,1,2 --trials 20000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s50_d12,ldpcdense_s50_d12_s2,ldpcdense_n12_s50_d12,ldpcdense_n12_s50_d12_s2,ldpcdense_n14_s50_d12,ldpcdense_s50_d15
run_csv "$RESULT_DIR/e5_K3200.csv" "$SIM" --K 3200 --oh 0,1,2 --trials 20000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s62_d12,ldpcdense_s62_d12_s2,ldpcdense_n12_s62_d12,ldpcdense_n12_s62_d12_s2,ldpcdense_n14_s62_d12,ldpcdense_s62_d19
run_csv "$RESULT_DIR/e5_K10000.csv" "$SIM" --K 10000 --oh 0,1,2 --trials 20000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s86_d12,ldpcdense_s86_d12_s2,ldpcdense_n12_s86_d12,ldpcdense_n12_s86_d12_s2,ldpcdense_n14_s86_d12,ldpcdense_s86_d26
mark_done /tmp/e5small_done
run_csv "$RESULT_DIR/e3a_K10000.csv" "$SIM" --K 10000 --oh 0,1,2,5 --trials 8000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s52_d12,ldpcdense_s52_d17,ldpcdense_s52_d26,ldpcdense_s60_d12,ldpcdense_s60_d17,ldpcdense_s60_d26,ldpcdense_s69_d12,ldpcdense_s69_d17,ldpcdense_s69_d26,ldpcdense_s77_d12,ldpcdense_s77_d17,ldpcdense_s77_d26,ldpcdense_s86_d12,ldpcdense_s86_d17,ldpcdense_s86_d26,ldpcdense_s69_d17_h12,ldpcdense_s86_d12_h12
run_csv "$RESULT_DIR/e3a_K16000.csv" "$SIM" --K 16000 --oh 0,1,2,5 --trials 8000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s68_d12,ldpcdense_s68_d23,ldpcdense_s68_d34,ldpcdense_s80_d12,ldpcdense_s80_d23,ldpcdense_s80_d34,ldpcdense_s91_d12,ldpcdense_s91_d23,ldpcdense_s91_d34,ldpcdense_s103_d12,ldpcdense_s103_d23,ldpcdense_s103_d34,ldpcdense_s114_d12,ldpcdense_s114_d23,ldpcdense_s114_d34,ldpcdense_s91_d23_h12,ldpcdense_s114_d12_h12
run_csv "$RESULT_DIR/e3a_K20000.csv" "$SIM" --K 20000 --oh 0,1,2,5 --trials 8000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s80_d12,ldpcdense_s80_d27,ldpcdense_s80_d40,ldpcdense_s94_d12,ldpcdense_s94_d27,ldpcdense_s94_d40,ldpcdense_s107_d12,ldpcdense_s107_d27,ldpcdense_s107_d40,ldpcdense_s121_d12,ldpcdense_s121_d27,ldpcdense_s121_d40,ldpcdense_s134_d12,ldpcdense_s134_d27,ldpcdense_s134_d40,ldpcdense_s107_d27_h12,ldpcdense_s134_d12_h12
run_csv "$RESULT_DIR/e3a_K32000.csv" "$SIM" --K 32000 --oh 0,1,2,5 --trials 5000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s114_d12,ldpcdense_s114_d38,ldpcdense_s114_d57,ldpcdense_s133_d12,ldpcdense_s133_d38,ldpcdense_s133_d57,ldpcdense_s152_d12,ldpcdense_s152_d38,ldpcdense_s152_d57,ldpcdense_s171_d12,ldpcdense_s171_d38,ldpcdense_s171_d57,ldpcdense_s190_d12,ldpcdense_s190_d38,ldpcdense_s190_d57,ldpcdense_s152_d38_h12,ldpcdense_s190_d12_h12
run_csv "$RESULT_DIR/e3a_K48000.csv" "$SIM" --K 48000 --oh 0,1,2,5 --trials 3000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s222_d12,ldpcdense_s222_d74,ldpcdense_s222_d111,ldpcdense_s259_d12,ldpcdense_s259_d74,ldpcdense_s259_d111,ldpcdense_s296_d12,ldpcdense_s296_d74,ldpcdense_s296_d111,ldpcdense_s333_d12,ldpcdense_s333_d74,ldpcdense_s333_d111,ldpcdense_s370_d12,ldpcdense_s370_d74,ldpcdense_s370_d111,ldpcdense_s296_d74_h12,ldpcdense_s370_d12_h12
run_csv "$RESULT_DIR/e3a_K64000.csv" "$SIM" --K 64000 --oh 0,1,2,5 --trials 3000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpc,ldpcdense_s208_d12,ldpcdense_s208_d69,ldpcdense_s208_d104,ldpcdense_s242_d12,ldpcdense_s242_d69,ldpcdense_s242_d104,ldpcdense_s277_d12,ldpcdense_s277_d69,ldpcdense_s277_d104,ldpcdense_s311_d12,ldpcdense_s311_d69,ldpcdense_s311_d104,ldpcdense_s346_d12,ldpcdense_s346_d69,ldpcdense_s346_d104,ldpcdense_s277_d69_h12,ldpcdense_s346_d12_h12
mark_done /tmp/e3abig_done
run_csv "$RESULT_DIR/e1_hgrid_K32000.csv" "$SIM" --K 32000 --oh 0,1,2 --trials 3000 --threads 116 \
  --paired --seed 24414209 --ge-replay \
  --schemes dense,dense_h12,ldpc,ldpc_h12,ldpcdense_s152_d38,ldpcdense_s152_d38_h8,ldpcdense_s152_d38_h12,ldpcdense_s190_d12_h12
run_csv "$RESULT_DIR/e1_hgrid_K64000.csv" "$SIM" --K 64000 --oh 0,1,2 --trials 1500 --threads 116 \
  --paired --seed 24414209 --ge-replay \
  --schemes dense,dense_h12,ldpc,ldpc_h12,ldpcdense_s277_d69,ldpcdense_s277_d69_h8,ldpcdense_s277_d69_h12,ldpcdense_s346_d12_h12
mark_done /tmp/e1bbig_done
run_csv "$RESULT_DIR/e5_K32000.csv" "$SIM" --K 32000 --oh 0,1,2 --trials 5000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpcdense_s190_d12,ldpcdense_s190_d12_s2,ldpcdense_n12_s190_d12,ldpcdense_n12_s190_d12_s2,ldpcdense_s190_d57
run_csv "$RESULT_DIR/e5_K64000.csv" "$SIM" --K 64000 --oh 0,1,2 --trials 3000 --threads 116 \
  --paired --seed 24414209 \
  --schemes dense,ldpcdense_s346_d12,ldpcdense_s346_d12_s2,ldpcdense_n12_s346_d12,ldpcdense_n12_s346_d12_s2,ldpcdense_s346_d104
mark_done /tmp/e5big_done
mark_done /tmp/queue_done
