#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/wirehair-test-verify-fix.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

cat >"$TMP/fake-whx" <<'EOF'
#!/bin/bash
set -euo pipefail
mode=${1:-}; shift || true
case "$mode" in
    scan)
        nfile=
        while [ "$#" -gt 0 ]; do
            if [ "$1" = --nfile ]; then nfile=$2; shift 2; else shift; fi
        done
        count=$(awk '/^[[:space:]]*(#|$)/ {next} {++n} END {print n+0}' "$nfile")
        if [ "${SCENARIO:-pass}" = malformed ]; then
            echo "not scan output"
            exit 0
        fi
        echo "# scan: fake"
        if [ "${SCENARIO:-pass}" = residual ]; then
            echo "WEAK N=100 mean=1.0000 max=1 fail=0"
        fi
        echo "# scan done: $count N scanned"
        ;;
    fuzz)
        echo "# fuzz: fake"
        if [ "${SCENARIO:-pass}" = fuzz_failure ]; then
            echo "# fuzz done: 0.0s trials=1 fails=1 rate=1 trials/s"
            echo "# overhead(success trials=0): mean=0 p50=0 p99=0 p999=0 max=0"
            exit 7
        fi
        echo "# fuzz done: 0.0s trials=1 fails=0 rate=1 trials/s"
        echo "# overhead(success trials=1): mean=0 p50=0 p99=0 p999=0 max=0"
        ;;
    ohead)
        value=0
        case "$(basename "$0"):${SCENARIO:-pass}" in
            *baseline*:*) value=0 ;;
            *:collateral) value=1 ;;
        esac
        echo "# ohead: fake"
        echo "N mean p50 p99 p999 max fail"
        echo "64 $value 0 0 0 0 0"
        ;;
    *) exit 2 ;;
esac
EOF
chmod +x "$TMP/fake-whx"
cp -f "$TMP/fake-whx" "$TMP/fake-baseline"

cat >"$TMP/fake-build" <<'EOF'
#!/bin/bash
set -euo pipefail
if [ "${SCENARIO:-pass}" = build_failure ]; then
    echo "injected build failure" >&2
    exit 9
fi
cp -f "$FAKE_HARNESS" "$OUT"
if [ "${SCENARIO:-pass}" != same_binary ]; then
    printf '\n# fixed-build marker\n' >>"$OUT"
fi
chmod +x "$OUT"
EOF
chmod +x "$TMP/fake-build"

printf '100\n200\n' >"$TMP/weak.txt"

run_case() {
    scenario=$1
    expected_rc=$2
    diagnostic=$3
    rm -f "$TMP/fixed"
    set +e
    SCENARIO="$scenario" FAKE_HARNESS="$TMP/fake-whx" \
        VERIFY_BUILD_SCRIPT="$TMP/fake-build" \
        VERIFY_FIXED_BIN="$TMP/fixed" \
        VERIFY_BASELINE_BIN="$TMP/fake-baseline" \
        VERIFY_WEAK_N_FILE="$TMP/weak.txt" \
        VERIFY_SCAN_TRIALS=1 VERIFY_FUZZ_SECONDS=0 VERIFY_FUZZ_NMAX=2 \
        VERIFY_COLLATERAL_TRIALS=1 VERIFY_COLLATERAL_NS=64 VERIFY_THREADS=1 \
        bash bench/verify_fix.sh >"$TMP/$scenario.out" 2>"$TMP/$scenario.err"
    rc=$?
    set -e
    if [ "$rc" -ne "$expected_rc" ]; then
        cat "$TMP/$scenario.out" "$TMP/$scenario.err" >&2
        echo "case $scenario: expected rc=$expected_rc, got rc=$rc" >&2
        exit 1
    fi
    if ! grep -Fq "$diagnostic" "$TMP/$scenario.out" "$TMP/$scenario.err"; then
        cat "$TMP/$scenario.out" "$TMP/$scenario.err" >&2
        echo "case $scenario: missing diagnostic: $diagnostic" >&2
        exit 1
    fi
}

run_case build_failure 1 "fixed harness build failed"
run_case same_binary 1 "baseline and fixed harness binaries are byte-identical"
run_case residual 1 "residual weak rows remain"
run_case fuzz_failure 1 "fuzz command failed"
run_case collateral 1 "non-fixed collateral output changed"
run_case malformed 1 "malformed residual scan output"
run_case pass 0 "verify_fix: PASS"

echo "test_verify_fix: PASS (7 scenarios)"
