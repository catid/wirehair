#!/bin/bash
# Assert the fixup contract against an exact weak set and a saved pre-fix binary.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

die() {
    echo "verify_fix: $*" >&2
    exit 1
}

THREADS=${VERIFY_THREADS:-4}
SCAN_TRIALS=${VERIFY_SCAN_TRIALS:-200}
SCAN_BB=${VERIFY_SCAN_BB:-64}
SCAN_STARTMODE=${VERIFY_SCAN_STARTMODE:-0}
SCAN_LOSS=${VERIFY_SCAN_LOSS:-0.10}
SCAN_THRESHOLD=${VERIFY_SCAN_THRESHOLD:-0.10}
SCAN_SEED=${VERIFY_SCAN_SEED:-0x5CA4}
FUZZ_SECONDS=${VERIFY_FUZZ_SECONDS:-30}
FUZZ_NMAX=${VERIFY_FUZZ_NMAX:-24000}
FUZZ_SEED=${VERIFY_FUZZ_SEED:-0xF1}
COLLATERAL_TRIALS=${VERIFY_COLLATERAL_TRIALS:-200}
COLLATERAL_SEED=${VERIFY_COLLATERAL_SEED:-0xC011A7}
COLLATERAL_NS=${VERIFY_COLLATERAL_NS:-64,128,512,1024,2048}
WEAK_N_FILE=${VERIFY_WEAK_N_FILE:-bench/results/weakN.txt}
FIXED_BIN=${VERIFY_FIXED_BIN:-bench/whx.fix}
BASELINE_BIN=${VERIFY_BASELINE_BIN:-bench/whx.baseline}
BUILD_SCRIPT=${VERIFY_BUILD_SCRIPT:-bench/build.sh}
SKIP_BUILD=${VERIFY_SKIP_BUILD:-0}

if [ -n "${VERIFY_WORKDIR:-}" ]; then
    WORK=$VERIFY_WORKDIR
    mkdir -p "$WORK"
    CLEAN_WORK=0
else
    WORK=$(mktemp -d "${TMPDIR:-/tmp}/wirehair-verify-fix.XXXXXX")
    CLEAN_WORK=1
fi
cleanup() {
    if [ "$CLEAN_WORK" -eq 1 ]; then rm -rf "$WORK"; fi
}
trap cleanup EXIT

case "$THREADS:$SCAN_TRIALS:$FUZZ_NMAX:$COLLATERAL_TRIALS" in
    *[!0-9:]*|0:*|*:0:*|*:*:*:0) die "thread/trial/range settings must be positive integers" ;;
esac
[ -r "$WEAK_N_FILE" ] || die "weak-set file is not readable: $WEAK_N_FILE"
[ -x "$BASELINE_BIN" ] || die "saved pre-fix baseline is not executable: $BASELINE_BIN"
[ "$BASELINE_BIN" != "$FIXED_BIN" ] || die "baseline and fixed harness paths must differ"
if [ -e "$FIXED_BIN" ] && [ "$BASELINE_BIN" -ef "$FIXED_BIN" ]; then
    die "baseline and fixed harness paths resolve to the same file"
fi

if [ "$SKIP_BUILD" != 1 ]; then
    echo "=== rebuild fixed harness ==="
    rm -f "$FIXED_BIN"
    if ! OUT="$FIXED_BIN" TAG=verify-fix bash "$BUILD_SCRIPT" >"$WORK/build.log" 2>&1; then
        cat "$WORK/build.log" >&2
        die "fixed harness build failed"
    fi
fi
[ -x "$FIXED_BIN" ] || die "fixed harness is not executable: $FIXED_BIN"
if [ "$BASELINE_BIN" -ef "$FIXED_BIN" ]; then
    die "build produced the baseline file instead of a distinct fixed harness"
fi
if cmp -s "$BASELINE_BIN" "$FIXED_BIN"; then
    die "baseline and fixed harness binaries are byte-identical"
fi

EXPECTED_COUNT=$(awk '
    /^[[:space:]]*(#|$)/ { next }
    $0 !~ /^[[:space:]]*[0-9]+[[:space:]]*$/ { exit 2 }
    { ++count }
    END { if (count > 0) print count; else exit 2 }
' "$WEAK_N_FILE") || die "weak-set file is empty or malformed: $WEAK_N_FILE"

echo "=== exact residual scan ($EXPECTED_COUNT requested N) ==="
if ! "$FIXED_BIN" scan --threads "$THREADS" --nfile "$WEAK_N_FILE" \
        --trials "$SCAN_TRIALS" --bb "$SCAN_BB" --startmode "$SCAN_STARTMODE" \
        --loss "$SCAN_LOSS" --thresh "$SCAN_THRESHOLD" --seed "$SCAN_SEED" \
        >"$WORK/residual.log" 2>"$WORK/residual.err"; then
    cat "$WORK/residual.err" >&2
    die "exact residual scan command failed"
fi

SCAN_COUNT=$(awk '
    /^# scan:/ { ++headers; next }
    /^WEAK N=/ { ++weak; next }
    /^# scan done: [0-9]+ N scanned$/ { ++summaries; scanned=$4; next }
    NF { bad=1 }
    END {
        if (headers != 1 || summaries != 1 || bad) exit 2
        print scanned
    }
' "$WORK/residual.log") || die "malformed residual scan output"
[ "$SCAN_COUNT" = "$EXPECTED_COUNT" ] || \
    die "residual scan covered $SCAN_COUNT N, expected exactly $EXPECTED_COUNT"
if grep -q '^WEAK N=' "$WORK/residual.log"; then
    grep '^WEAK N=' "$WORK/residual.log" >&2
    die "residual weak rows remain"
fi

echo "=== exact recovery fuzz ==="
if ! "$FIXED_BIN" fuzz --threads "$THREADS" --secs "$FUZZ_SECONDS" \
        --nmax "$FUZZ_NMAX" --seed "$FUZZ_SEED" \
        >"$WORK/fuzz.log" 2>"$WORK/fuzz.err"; then
    cat "$WORK/fuzz.err" >&2
    tail -20 "$WORK/fuzz.log" >&2 || true
    die "fuzz command failed"
fi
awk '
    /^# fuzz:/ { ++headers; next }
    /^HIGH-OVERHEAD / { next }
    /^# fuzz done: .* fails=0 rate=/ { ++summaries; next }
    /^# overhead\(success trials=[0-9]+\):/ { ++overheads; next }
    NF { bad=1 }
    END { if (headers != 1 || summaries != 1 || overheads != 1 || bad) exit 2 }
' "$WORK/fuzz.log" || die "malformed fuzz output or nonzero reported failures"

COLLATERAL_FILE="$WORK/collateral-n.txt"
awk -v list="$COLLATERAL_NS" 'BEGIN {
    count = split(list, values, ",")
    for (i = 1; i <= count; ++i) {
        if (values[i] !~ /^[0-9]+$/ || values[i] + 0 < 2 || values[i] + 0 > 64000) exit 2
        print values[i] + 0
    }
}' >"$COLLATERAL_FILE" || die "VERIFY_COLLATERAL_NS is malformed"

while IFS= read -r N; do
    if grep -Eq "\{[[:space:]]*${N}[[:space:]]*," WirehairPeelFixups.inc WirehairDenseFixups.inc; then
        die "collateral N=$N is present in a shipped fixup table"
    fi
done <"$COLLATERAL_FILE"

echo "=== deterministic non-fixed collateral comparison ==="
: >"$WORK/collateral.baseline"
: >"$WORK/collateral.fixed"
while IFS= read -r N; do
    for SPEC in "baseline:$BASELINE_BIN" "fixed:$FIXED_BIN"; do
        LABEL=${SPEC%%:*}
        BIN=${SPEC#*:}
        if ! "$BIN" ohead --threads "$THREADS" --nlo "$N" --nhi "$N" --nstep 1 \
                --trials "$COLLATERAL_TRIALS" --bb "$SCAN_BB" \
                --startmode "$SCAN_STARTMODE" --loss "$SCAN_LOSS" \
                --seed "$COLLATERAL_SEED" >>"$WORK/collateral.$LABEL" 2>"$WORK/collateral.$LABEL.err"; then
            cat "$WORK/collateral.$LABEL.err" >&2
            die "$LABEL collateral command failed for N=$N"
        fi
    done
done <"$COLLATERAL_FILE"

COLLATERAL_COUNT=$(wc -l <"$COLLATERAL_FILE")
for LABEL in baseline fixed; do
    awk -v expected="$COLLATERAL_COUNT" '
        /^# ohead:/ { ++headers; next }
        /^N[[:space:]]+mean[[:space:]]+p50[[:space:]]+p99[[:space:]]+p999[[:space:]]+max[[:space:]]+fail[[:space:]]*$/ { ++columns; next }
        $1 ~ /^[0-9]+$/ && NF == 7 { ++rows; next }
        NF { bad=1 }
        END { if (headers != expected || columns != expected || rows != expected || bad) exit 2 }
    ' "$WORK/collateral.$LABEL" || die "malformed $LABEL collateral output"
done

if ! cmp -s "$WORK/collateral.baseline" "$WORK/collateral.fixed"; then
    diff -u "$WORK/collateral.baseline" "$WORK/collateral.fixed" >&2 || true
    die "non-fixed collateral output changed"
fi

echo "verify_fix: PASS exact_N=$EXPECTED_COUNT fuzz_seed=$FUZZ_SEED collateral=$COLLATERAL_COUNT"
