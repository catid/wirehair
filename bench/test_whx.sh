#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/wirehair-test-whx.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

NORMAL="$TMP/whx"
KNOBS="$TMP/whx-knobs"
OUT="$NORMAL" TAG=test-whx-normal bash bench/build.sh >/dev/null
OUT="$KNOBS" TAG=test-whx-knobs EXTRA=-DWH_SEED_KNOBS=1 bash bench/build.sh >/dev/null

"$NORMAL" selftest
"$KNOBS" selftest

"$NORMAL" ohead --threads 1 --nlo 100 --nhi 100 --trials 1 --bb 1 \
    --startmode 1 --loss 0 --seed 12 >"$TMP/ohead-one.out"
awk '$1 == 100 { found=1; good=($2 == 1 && $3 == 1 && $4 == 1 && $5 == 1 && $6 == 1 && $7 == 0) }
     END { exit !(found && good) }' "$TMP/ohead-one.out"

"$NORMAL" ohead --threads 3 --nlo 96 --nhi 96 --trials 41 --bb 7 \
    --startmode 0 --loss 0.17 --seed 12345 --samples-out "$TMP/ohead.samples" \
    >"$TMP/ohead.out"
reference_quantile() {
    probability=$1
    count=$(awk '!/^#/ {++n} END {print n+0}' "$TMP/ohead.samples")
    rank=$(awk -v p="$probability" -v n="$count" 'BEGIN { x=p*n; r=int(x); if (r<x) ++r; if (r<1) r=1; print r }')
    awk '!/^#/ {print $3}' "$TMP/ohead.samples" | sort -n | sed -n "${rank}p"
}
read -r _ _ reported_p50 reported_p99 reported_p999 _ reported_fail < <(awk '$1 == 96 {print}' "$TMP/ohead.out")
[ "$reported_fail" = 0 ]
[ "$reported_p50" = "$(reference_quantile 0.50)" ]
[ "$reported_p99" = "$(reference_quantile 0.99)" ]
[ "$reported_p999" = "$(reference_quantile 0.999)" ]
"$NORMAL" ohead --threads 1 --nlo 32 --nhi 64 --nstep 2147483647 --trials 1 \
    --bb 1 --startmode 0 --loss 0 --seed 9 >"$TMP/ohead-large-step.out"
grep -Eq '^32[[:space:]]' "$TMP/ohead-large-step.out"

"$NORMAL" fuzz --threads 2 --secs 0.2 --nmax 128 --seed 777 >"$TMP/fuzz.out"
grep -Eq '^# fuzz done: .* fails=0 rate=' "$TMP/fuzz.out"

"$NORMAL" bench --threads 1 --N 4096 --bb 64 --rounds 1 --memory-mib 64 >"$TMP/bench-one.out"
"$NORMAL" bench --threads 64 --N 4096 --bb 64 --rounds 1 --memory-mib 64 >"$TMP/bench-many.out"
grep -Eq 'workers=1 .*verified=1' "$TMP/bench-one.out"
grep -Eq 'workers=1 .*verified=1' "$TMP/bench-many.out"
one_overhead=$(awk '$1 == 4096 {print $6}' "$TMP/bench-one.out")
many_overhead=$(awk '$1 == 4096 {print $6}' "$TMP/bench-many.out")
[ "$one_overhead" = "$many_overhead" ]

"$NORMAL" bench --threads 4 --N 512 --bb 64 --rounds 4 --memory-mib 64 >"$TMP/bench-scale.out"
grep -Eq 'workers=4 .*verified=4' "$TMP/bench-scale.out"

set +e
"$NORMAL" bench --threads 4 --N 4096 --bb 64 --rounds 2 --memory-mib 0 \
    >"$TMP/memory.out" 2>"$TMP/memory.err"
memory_rc=$?
set -e
[ "$memory_rc" -eq 2 ]
grep -Fq 'positive --memory-mib' "$TMP/memory.err"

for location in first middle last; do
    set +e
    WHX_BENCH_CORRUPT=$location "$NORMAL" bench --threads 4 --N 32 --bb 7 \
        --rounds 1 --memory-mib 1 >"$TMP/fault-$location.out" 2>"$TMP/fault-$location.err"
    fault_rc=$?
    set -e
    [ "$fault_rc" -eq 1 ]
    grep -Fq 'BENCH_MISMATCH seed=' "$TMP/fault-$location.err"
    grep -Eq 'FAILS=1 MISMATCHES=1 .*verified=0' "$TMP/fault-$location.out"
done
grep -Fq 'offset=0 ' "$TMP/fault-first.err"
grep -Fq 'offset=112 ' "$TMP/fault-middle.err"
grep -Fq 'offset=223 ' "$TMP/fault-last.err"
set +e
WHX_BENCH_CORRUPT=999 "$NORMAL" bench --threads 1 --N 32 --bb 7 --rounds 1 \
    --memory-mib 1 >"$TMP/fault-range.out" 2>"$TMP/fault-range.err"
fault_range_rc=$?
set -e
[ "$fault_range_rc" -eq 2 ]
grep -Fq 'is outside N=32 bb=7 message' "$TMP/fault-range.err"

printf '32\n64\n' >"$TMP/n-file"
"$NORMAL" scan --threads 64 --nfile "$TMP/n-file" --trials 1 --bb 1 \
    --startmode 0 --loss 0.10 --thresh 999 --seed 1 >"$TMP/scan.out"
grep -Fq '# scan done: 2 N scanned' "$TMP/scan.out"

for bad_loss in 0.20 nan malformed; do
    set +e
    "$KNOBS" seedsearch --threads 1 --nlist 32 --nseeds 1 --tsearch 1 --tverify 1 \
        --loss "$bad_loss" >"$TMP/loss-$bad_loss.out" 2>"$TMP/loss-$bad_loss.err"
    loss_rc=$?
    set -e
    [ "$loss_rc" -eq 2 ]
    [ ! -s "$TMP/loss-$bad_loss.out" ]
done

seedsearch_args=(seedsearch --nlist 32 --nseeds 3 --dseeds 2 --goodthr 0
                 --tsearch 2 --tverify 4 --bb 0 --loss 0.10)
"$KNOBS" "${seedsearch_args[@]}" --threads 1 >"$TMP/seed-t1.out" 2>"$TMP/seed-t1.err"
"$KNOBS" "${seedsearch_args[@]}" --threads 4 >"$TMP/seed-t4.out" 2>"$TMP/seed-t4.err"
cmp "$TMP/seed-t1.out" "$TMP/seed-t4.out"
grep -Fq 'primary_loss=0.10 secondary_loss=0.30 finalist_trials=paired' "$TMP/seed-t1.out"
grep -Fq 'pairing=paired validation_trials=4 trial_stride=2654435761 primary_base=0xBEE5+N*577 secondary_base=0xF00D+N*331' "$TMP/seed-t1.out"

OUT="$TMP/peel.inc" DEFTHRESH=-1 GOOD=999 GOOD30=999 \
    bash bench/gen_fixups.sh "$TMP/seed-t1.out" >/dev/null
PEEL_OUT="$TMP/peel2.inc" DENSE_OUT="$TMP/dense2.inc" \
    DEFTHRESH=-1 GOOD=999 GOOD30=999 \
    bash bench/gen_fixups2.sh "$TMP/seed-t1.out" >/dev/null
cat >"$TMP/include-test.cpp" <<EOF
struct Fixup { int n; int seed; };
static const Fixup peel[] = {
#include "$TMP/peel.inc"
};
static const Fixup peel2[] = {
#include "$TMP/peel2.inc"
};
static const Fixup dense2[] = {
#include "$TMP/dense2.inc"
};
int main() { return peel[0].n < 1 || peel2[0].n < 1 || dense2[0].n < 1; }
EOF
g++ -std=c++11 -Wall -Wextra -Werror "$TMP/include-test.cpp" -o "$TMP/include-test"
"$TMP/include-test"

# Link the generated tables into the real codec by compiling a private copy of
# WirehairTools.cpp (quoted fixup includes resolve beside that copy).
cp -f WirehairTools.cpp "$TMP/WirehairTools.cpp"
CODEC_FLAGS=(-O2 -march=native -std=c++11 -pthread -Wall -Wno-unused-function
             -I"$ROOT" -I"$ROOT/include" -I"$ROOT/test")
g++ "${CODEC_FLAGS[@]}" -c "$TMP/WirehairTools.cpp" -o "$TMP/WirehairTools.generated.o"
g++ bench/obj/test-whx-normal/gf256.o \
    bench/obj/test-whx-normal/WirehairCodec.o "$TMP/WirehairTools.generated.o" \
    bench/obj/test-whx-normal/wirehair.o bench/obj/test-whx-normal/whx.o \
    -pthread -o "$TMP/whx-generated"
"$TMP/whx-generated" fuzz --threads 2 --secs 0.1 --nmax 64 --seed 8080 \
    >"$TMP/generated-fuzz.out"
grep -Eq '^# fuzz done: .* fails=0 rate=' "$TMP/generated-fuzz.out"
for generated_loss in 0.10 0.30; do
    "$TMP/whx-generated" ohead --threads 2 --nlo 32 --nhi 32 --trials 8 --bb 7 \
        --startmode 0 --loss "$generated_loss" --seed 8181 \
        >"$TMP/generated-ohead-$generated_loss.out"
    awk '$1 == 32 { found=1; good=($7 == 0) } END { exit !(found && good) }' \
        "$TMP/generated-ohead-$generated_loss.out"
done

read -r selected_n selected_p selected_d _ _ _ < <(awk '!/^#/ {print; exit}' "$TMP/seed-t1.out")
for loss_seed in '0.10 9001' '0.30 9002'; do
    read -r cert_loss cert_seed <<<"$loss_seed"
    "$KNOBS" seedmean --threads 2 --N "$selected_n" --pseed "$selected_p" --dseed "$selected_d" \
        --trials 8 --bb 0 --loss "$cert_loss" --startmode 0 --seed "$cert_seed" \
        >"$TMP/cert-$cert_loss.out"
    awk '$1 ~ /^[0-9]+$/ { found=1; good=($6 == 0) } END { exit !(found && good) }' \
        "$TMP/cert-$cert_loss.out"
done

printf '# N best_pseed default_mean best_mean@0.20 best_mean@0.30\n32 1 1 0 0\n' \
    >"$TMP/mislabeled.out"
printf 'sentinel\n' >"$TMP/preserved.inc"
set +e
OUT="$TMP/preserved.inc" bash bench/gen_fixups.sh "$TMP/mislabeled.out" \
    >"$TMP/mislabeled.stdout" 2>"$TMP/mislabeled.stderr"
generator_rc=$?
set -e
[ "$generator_rc" -ne 0 ]
[ "$(cat "$TMP/preserved.inc")" = sentinel ]

cat >"$TMP/mixed-width.out" <<'EOF'
# schema=whx-seedsearch-v1 primary_loss=0.10 secondary_loss=0.30 finalist_trials=paired tie_break=peel_then_lowest_seed
# pairing=paired validation_trials=4 trial_stride=2654435761 primary_base=0xBEE5+N*577 secondary_base=0xF00D+N*331
# N best_pseed best_dseed default_mean best_mean@0.10 best_mean@0.30
32 1 -1 0.1 0.1 0.1
33 1 0.1 0.1 0.1
EOF
printf 'sentinel\n' >"$TMP/preserved-mixed.inc"
set +e
OUT="$TMP/preserved-mixed.inc" bash bench/gen_fixups.sh "$TMP/mixed-width.out" \
    >"$TMP/mixed.stdout" 2>"$TMP/mixed.stderr"
mixed_rc=$?
set -e
[ "$mixed_rc" -ne 0 ]
[ "$(cat "$TMP/preserved-mixed.inc")" = sentinel ]
grep -Fq 'data column count does not match header' "$TMP/mixed.stderr"

bash bench/test_verify_fix.sh

if [ "${RUN_SANITIZERS:-0}" = 1 ]; then
    SAN="$TMP/whx-san"
    SAN_FLAGS='-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer'
    OUT="$SAN" TAG=test-whx-san EXTRA="$SAN_FLAGS" LINK_EXTRA='-fsanitize=address,undefined' \
        bash bench/build.sh >/dev/null
    ASAN_OPTIONS=detect_leaks=1 "$SAN" selftest
    ASAN_OPTIONS=detect_leaks=1 "$SAN" ohead --threads 2 --nlo 32 --nhi 32 \
        --trials 4 --bb 7 --startmode 0 --loss 0.10 --seed 55 >/dev/null
    ASAN_OPTIONS=detect_leaks=1 "$SAN" bench --threads 8 --N 32 --bb 7 \
        --rounds 2 --memory-mib 2 >/dev/null
fi

echo "test_whx: PASS"
