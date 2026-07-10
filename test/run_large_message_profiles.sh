#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "usage: $0 <large_message_test> [quick|scheduled|all]" >&2
    exit 2
}

[[ $# -ge 1 && $# -le 2 ]] || usage
binary=$1
mode=${2:-quick}
[[ -x "$binary" ]] || {
    echo "large-message runner: not executable: $binary" >&2
    exit 2
}
case "$mode" in
    quick|scheduled|all) ;;
    *) usage ;;
esac

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/wirehair-large-message.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT
have_gnu_time=0
if [[ -x /usr/bin/time ]] &&
    /usr/bin/time -f '' -o /dev/null true >/dev/null 2>&1
then
    have_gnu_time=1
fi

run_profile() {
    local name=$1
    local expectation=$2
    local expected_rc=$3
    local pattern=$4
    shift 4

    local output="$tmpdir/$name.log"
    local timing="$tmpdir/$name.time"
    local rc
    echo "== large-message profile: $name ($expectation) =="

    set +e
    if [[ $have_gnu_time -eq 1 ]]; then
        /usr/bin/time \
            -f "profile=$name elapsed_s=%e max_rss_kib=%M" \
            -o "$timing" \
            env "$@" "$binary" >"$output" 2>&1
        rc=$?
    else
        env "$@" "$binary" >"$output" 2>&1
        rc=$?
        printf 'profile=%s external_time_unavailable\n' "$name" >"$timing"
    fi
    set -e

    cat "$output"
    cat "$timing"
    if [[ "$expectation" == success ]]; then
        if [[ $rc -ne 0 ]] || ! grep -q 'large_message_test: PASS' "$output"; then
            echo "large-message runner: $name failed with status $rc" >&2
            return 1
        fi
    else
        if [[ $rc -ne $expected_rc ]] || ! grep -q "$pattern" "$output" ||
            grep -Eq 'AddressSanitizer|LeakSanitizer|UndefinedBehaviorSanitizer|MemorySanitizer|ThreadSanitizer|runtime error:' "$output"
        then
            echo "large-message runner: $name did not produce expected failure" >&2
            return 1
        fi
    fi
}

run_quick() {
    local common=(
        WIREHAIR_LARGE_N=256
        WIREHAIR_LARGE_BLOCK_BYTES=257
        WIREHAIR_LARGE_FINAL_BYTES=113
        WIREHAIR_LARGE_LOSS_COUNT=5
        WIREHAIR_LARGE_MAX_REPAIR=64
        WIREHAIR_LARGE_SEED=0x8af09d31e7c4526b
        WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=1048576
        WIREHAIR_LARGE_MAX_MILLISECONDS=30000
    )
    run_profile fast_partial success 0 '' "${common[@]}"
    run_profile fast_partial_repeat success 0 '' "${common[@]}"
    local first_digest_count
    local repeat_digest_count
    local first_digest
    local repeat_digest
    first_digest_count=$(grep -c '^packet_order_digest=.' \
        "$tmpdir/fast_partial.log" || true)
    repeat_digest_count=$(grep -c '^packet_order_digest=.' \
        "$tmpdir/fast_partial_repeat.log" || true)
    if [[ $first_digest_count -ne 1 || $repeat_digest_count -ne 1 ]]
    then
        echo "large-message runner: expected one nonempty packet digest per replay" >&2
        return 1
    fi
    first_digest=$(grep '^packet_order_digest=.' "$tmpdir/fast_partial.log")
    repeat_digest=$(grep '^packet_order_digest=.' \
        "$tmpdir/fast_partial_repeat.log")
    if [[ $first_digest != "$repeat_digest" ]]
    then
        echo "large-message runner: seeded packet order is not deterministic" >&2
        return 1
    fi
    run_profile terminal_decode failure 1 'Injected terminal decode result' \
        "${common[@]}" WIREHAIR_LARGE_INJECT_TERMINAL=1
    run_profile repair_cap failure 1 'repair blocks (cap reached)' \
        "${common[@]}" WIREHAIR_LARGE_LOSS_COUNT=8 \
        WIREHAIR_LARGE_MAX_REPAIR=1
    run_profile timeout_path failure 3 'Injected timeout failure' \
        "${common[@]}" WIREHAIR_LARGE_INJECT_TIMEOUT=1
    run_profile resource_cap failure 2 'Resource policy rejected' \
        "${common[@]}" WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=1024
}

run_scheduled() {
    run_profile n64000_byte success 0 '' \
        WIREHAIR_LARGE_N=64000 \
        WIREHAIR_LARGE_BLOCK_BYTES=1 \
        WIREHAIR_LARGE_FINAL_BYTES=1 \
        WIREHAIR_LARGE_LOSS_COUNT=17 \
        WIREHAIR_LARGE_MAX_REPAIR=128 \
        WIREHAIR_LARGE_SEED=0x640001 \
        WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=65536 \
        WIREHAIR_LARGE_MAX_MILLISECONDS=300000

    run_profile packet_partial success 0 '' \
        WIREHAIR_LARGE_N=4096 \
        WIREHAIR_LARGE_BLOCK_BYTES=1280 \
        WIREHAIR_LARGE_FINAL_BYTES=733 \
        WIREHAIR_LARGE_LOSS_COUNT=11 \
        WIREHAIR_LARGE_MAX_REPAIR=96 \
        WIREHAIR_LARGE_SEED=0x40961280 \
        WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=8388608 \
        WIREHAIR_LARGE_MAX_MILLISECONDS=180000

    run_profile large_block_partial success 0 '' \
        WIREHAIR_LARGE_N=32 \
        WIREHAIR_LARGE_BLOCK_BYTES=1048576 \
        WIREHAIR_LARGE_FINAL_BYTES=1048559 \
        WIREHAIR_LARGE_LOSS_COUNT=3 \
        WIREHAIR_LARGE_MAX_REPAIR=64 \
        WIREHAIR_LARGE_SEED=0x321048576 \
        WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=33554432 \
        WIREHAIR_LARGE_MAX_MILLISECONDS=180000
}

if [[ "$mode" == quick || "$mode" == all ]]; then
    run_quick
fi
if [[ "$mode" == scheduled || "$mode" == all ]]; then
    run_scheduled
fi

echo "large-message profile suite: PASS ($mode)"
