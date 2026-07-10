#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
runner="$repo_root/test/run_large_message_profiles.sh"
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/wirehair-large-runner-test.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

fake="$tmpdir/fake_large_message_test"
cat >"$fake" <<'FAKE'
#!/usr/bin/env bash
set -eu

profile=success
expected_rc=0
diagnostic='large_message_test: PASS'
if [[ ${WIREHAIR_LARGE_INJECT_TERMINAL:-0} == 1 ]]; then
    profile=terminal
    expected_rc=1
    diagnostic='Injected terminal decode result'
elif [[ ${WIREHAIR_LARGE_INJECT_TIMEOUT:-0} == 1 ]]; then
    profile=timeout
    expected_rc=3
    diagnostic='Injected timeout failure'
elif [[ ${WIREHAIR_LARGE_MAX_PAYLOAD_BYTES:-} == 1024 ]]; then
    profile=resource
    expected_rc=2
    diagnostic='Resource policy rejected'
elif [[ ${WIREHAIR_LARGE_MAX_REPAIR:-} == 1 ]]; then
    profile=repair
    expected_rc=1
    diagnostic='Decoder still needs more after 1 repair blocks (cap reached)'
fi

if [[ ${WIREHAIR_FAKE_OMIT_DIGEST:-0} != 1 ]]; then
    echo 'packet_order_digest=123456789 dropped_ids=1,2,3'
fi
echo "$diagnostic"
if [[ ${WIREHAIR_FAKE_SANITIZER_CASE:-} == "$profile" ]]; then
    echo 'ERROR: AddressSanitizer: injected self-test signature' >&2
fi
if [[ ${WIREHAIR_FAKE_WRONG_CASE:-} == "$profile" ]]; then
    exit 7
fi
exit "$expected_rc"
FAKE
chmod +x "$fake"

"$runner" "$fake" quick >"$tmpdir/pass.log" 2>&1
grep -q 'large-message profile suite: PASS (quick)' "$tmpdir/pass.log"

for profile in terminal repair timeout resource; do
    log="$tmpdir/wrong-$profile.log"
    if WIREHAIR_FAKE_WRONG_CASE=$profile "$runner" "$fake" quick >"$log" 2>&1; then
        echo "runner accepted the wrong exit status for $profile" >&2
        exit 1
    fi
    grep -q 'did not produce expected failure' "$log"
done

sanitizer_log="$tmpdir/sanitizer.log"
if WIREHAIR_FAKE_SANITIZER_CASE=terminal \
    "$runner" "$fake" quick >"$sanitizer_log" 2>&1
then
    echo 'runner accepted a sanitizer failure on an expected-failure path' >&2
    exit 1
fi
grep -q 'did not produce expected failure' "$sanitizer_log"

missing_digest_log="$tmpdir/missing-digest.log"
if WIREHAIR_FAKE_OMIT_DIGEST=1 \
    "$runner" "$fake" quick >"$missing_digest_log" 2>&1
then
    echo 'runner accepted successful replays with no packet digest' >&2
    exit 1
fi
grep -q 'expected one nonempty packet digest per replay' "$missing_digest_log"

echo 'large-message runner self-test: PASS'
