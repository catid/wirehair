#!/usr/bin/env bash
# Launch the sealed normalized-H15-v5 recovery holdout at 120x1 concurrency.
set -euo pipefail
export LC_ALL=C TZ=UTC PYTHONDONTWRITEBYTECODE=1

leader_wait_rc=
last_group_forced_kill=0
last_group_quiescent=1
last_process_forced_kill=0
last_process_quiescent=1

process_group_alive() {
    local pgid=$1
    [[ $pgid =~ ^[1-9][0-9]*$ ]] && kill -0 -- "-$pgid" 2>/dev/null
}

process_group_has_live_members() {
    local pgid=$1
    [[ $pgid =~ ^[1-9][0-9]*$ ]] || return 1
    ps -e -o pgid=,stat= | awk -v pgid="$pgid" \
        '$1 == pgid && $2 !~ /^Z/ { found = 1 } END { exit found ? 0 : 1 }'
}

process_is_live() {
    local pid=$1
    local state
    [[ $pid =~ ^[1-9][0-9]*$ ]] || return 1
    state=$(ps -o stat= -p "$pid" 2>/dev/null | tr -d '[:space:]')
    [[ -n $state && $state != Z* ]]
}

terminate_process_group() {
    local pgid=$1
    local checks=${2:-50}
    local delay=${3:-0.1}
    local iteration
    last_group_forced_kill=0
    last_group_quiescent=1
    if ! process_group_alive "$pgid"; then
        return 0
    fi
    kill -TERM -- "-$pgid" 2>/dev/null || true
    for ((iteration = 0; iteration < checks; ++iteration)); do
        if ! process_group_has_live_members "$pgid"; then
            return 0
        fi
        sleep "$delay"
    done
    if process_group_alive "$pgid"; then
        kill -KILL -- "-$pgid" 2>/dev/null || true
        last_group_forced_kill=1
    fi
    for ((iteration = 0; iteration < checks; ++iteration)); do
        if ! process_group_has_live_members "$pgid"; then
            return 0
        fi
        sleep "$delay"
    done
    if process_group_has_live_members "$pgid"; then
        last_group_quiescent=0
    fi
}

terminate_process() {
    local pid=$1
    local checks=${2:-50}
    local delay=${3:-0.1}
    local iteration
    last_process_forced_kill=0
    last_process_quiescent=1
    if ! process_is_live "$pid"; then
        return 0
    fi
    kill -TERM "$pid" 2>/dev/null || true
    for ((iteration = 0; iteration < checks; ++iteration)); do
        if ! process_is_live "$pid"; then
            return 0
        fi
        sleep "$delay"
    done
    if process_is_live "$pid"; then
        kill -KILL "$pid" 2>/dev/null || true
        last_process_forced_kill=1
    fi
    for ((iteration = 0; iteration < checks; ++iteration)); do
        if ! process_is_live "$pid"; then
            return 0
        fi
        sleep "$delay"
    done
    if process_is_live "$pid"; then
        last_process_quiescent=0
    fi
}

wait_for_leader() {
    local pid=$1
    if wait "$pid" 2>/dev/null; then
        leader_wait_rc=0
    else
        leader_wait_rc=$?
    fi
}

cleanup() {
    local rc=$?
    local pgid=${xargs_pgid:-${xargs_pid:-}}
    local checks=${cleanup_checks:-50}
    local delay=${cleanup_delay:-0.1}
    trap '' INT TERM
    trap - EXIT
    if [[ -n $pgid ]]; then
        terminate_process_group "$pgid" "$checks" "$delay" || true
    fi
    if [[ -n ${xargs_pid:-} ]]; then
        wait_for_leader "$xargs_pid"
        xargs_pid=
    fi
    if [[ -n ${watchdog_pid:-} ]]; then
        terminate_process "$watchdog_pid" "$checks" "$delay" || true
        wait_for_leader "$watchdog_pid"
        watchdog_pid=
    fi
    if [[ -n ${launch_tmp:-} ]]; then
        rm -rf -- "$launch_tmp"
    fi
    return "$rc"
}

launcher_selftest_cleanup() {
    if [[ -n ${selftest_pgid:-} ]]; then
        terminate_process_group "$selftest_pgid" 2 0.05 || true
    fi
    if [[ -n ${selftest_tmp:-} ]]; then
        rm -rf -- "$selftest_tmp"
    fi
}

launcher_selftest() {
    local cleanup_rc leader stubborn_pid state
    selftest_tmp=$(mktemp -d -- /tmp/wh2-launch-selftest.XXXXXX)
    selftest_pgid=
    trap launcher_selftest_cleanup EXIT

    setsid bash -c 'exit 0' &
    leader=$!
    wait_for_leader "$leader"
    if [[ $leader_wait_rc != 0 ]]; then
        echo "launcher selftest lost normal leader rc" >&2
        return 1
    fi
    setsid bash -c 'exit 23' &
    leader=$!
    wait_for_leader "$leader"
    if [[ $leader_wait_rc != 23 ]]; then
        echo "launcher selftest lost failing leader rc" >&2
        return 1
    fi

    # shellcheck disable=SC2016
    setsid bash -c '
        trap "exit 0" TERM
        printf "ready\n" >"$1"
        while :; do sleep 1; done
    ' _ "$selftest_tmp/cooperative.ready" >/dev/null 2>&1 &
    leader=$!
    selftest_pgid=$leader
    while [[ ! -s $selftest_tmp/cooperative.ready ]]; do sleep 0.01; done
    terminate_process_group "$selftest_pgid" 5 0.05
    if [[ $last_group_forced_kill != 0 || $last_group_quiescent != 1 ]]; then
        echo "launcher selftest failed cooperative TERM group cleanup" >&2
        return 1
    fi
    wait_for_leader "$leader"
    selftest_pgid=
    if [[ $leader_wait_rc != 0 ]]; then
        echo "launcher selftest changed cooperative TERM leader rc" >&2
        return 1
    fi

    python3 -c '
import pathlib, signal, sys, time
signal.signal(signal.SIGTERM, signal.SIG_IGN)
pathlib.Path(sys.argv[1]).write_text("ready\n", encoding="ascii")
while True: time.sleep(1)
' "$selftest_tmp/direct.ready" &
    leader=$!
    while [[ ! -s $selftest_tmp/direct.ready ]]; do sleep 0.01; done
    terminate_process "$leader" 3 0.05
    if [[ $last_process_forced_kill != 1 || $last_process_quiescent != 1 ]]; then
        echo "launcher selftest failed direct TERM/KILL cleanup" >&2
        return 1
    fi
    wait_for_leader "$leader"
    if [[ $leader_wait_rc != 137 ]]; then
        echo "launcher selftest changed SIGKILLed direct leader rc" >&2
        return 1
    fi

    # Positional parameters intentionally expand in these isolated shells.
    # shellcheck disable=SC2016
    setsid bash -c '
        bash -c '\''trap "" TERM; printf "ready\\n" >"$1"; while :; do sleep 1; done'\'' \
            _ "$2" </dev/null >/dev/null 2>&1 &
        printf "%s\n" "$!" >"$1"
        while [[ ! -s $2 ]]; do sleep 0.01; done
        exit 7
    ' _ "$selftest_tmp/stubborn.pid" "$selftest_tmp/stubborn.ready" &
    leader=$!
    selftest_pgid=$leader
    wait_for_leader "$leader"
    if [[ $leader_wait_rc != 7 || ! -s $selftest_tmp/stubborn.pid ]]; then
        echo "launcher selftest did not reproduce exited leader topology" >&2
        return 1
    fi
    stubborn_pid=$(<"$selftest_tmp/stubborn.pid")
    if ! process_group_alive "$selftest_pgid" || ! kill -0 "$stubborn_pid" 2>/dev/null; then
        echo "launcher selftest stubborn descendant was not alive after leader exit" >&2
        return 1
    fi
    terminate_process_group "$selftest_pgid" 3 0.05
    if [[ $last_group_forced_kill != 1 || $last_group_quiescent != 1 ]]; then
        echo "launcher selftest did not KILL the SIGTERM-ignoring group" >&2
        return 1
    fi
    while read -r state; do
        if [[ -n $state && $state != Z* ]]; then
            echo "launcher selftest left a live process-group member: $state" >&2
            return 1
        fi
    done < <(ps -e -o pgid=,stat= | awk -v pgid="$selftest_pgid" '$1 == pgid { print $2 }')
    selftest_pgid=

    if (
        xargs_pid=
        xargs_pgid=
        watchdog_pid=
        launch_tmp=$selftest_tmp/signal-cleanup-tmp
        cleanup_checks=3
        cleanup_delay=0.05
        mkdir -m 700 -- "$launch_tmp"
        setsid python3 -c '
import pathlib, signal, sys, time
signal.signal(signal.SIGTERM, signal.SIG_IGN)
pathlib.Path(sys.argv[1]).write_text("ready\n", encoding="ascii")
while True: time.sleep(1)
' "$selftest_tmp/cleanup-xargs.ready" &
        xargs_pid=$!
        xargs_pgid=$xargs_pid
        python3 -c '
import pathlib, signal, sys, time
signal.signal(signal.SIGTERM, signal.SIG_IGN)
pathlib.Path(sys.argv[1]).write_text("ready\n", encoding="ascii")
while True: time.sleep(1)
' "$selftest_tmp/cleanup-watchdog.ready" &
        watchdog_pid=$!
        printf '%s\n' "$xargs_pid" >"$selftest_tmp/cleanup-xargs.pid"
        printf '%s\n' "$watchdog_pid" >"$selftest_tmp/cleanup-watchdog.pid"
        while [[ ! -s $selftest_tmp/cleanup-xargs.ready ||
                 ! -s $selftest_tmp/cleanup-watchdog.ready ]]; do sleep 0.01; done
        trap 'cleanup 2>/dev/null' EXIT
        trap 'exit 143' TERM
        target=$BASHPID
        (
            sleep 0.02
            kill -TERM "$target"
            sleep 0.03
            kill -TERM "$target" 2>/dev/null || true
        ) &
        while :; do sleep 1; done
    ); then
        cleanup_rc=0
    else
        cleanup_rc=$?
    fi
    if [[ $cleanup_rc != 143 || -e $selftest_tmp/signal-cleanup-tmp ]]; then
        echo "launcher selftest signal-driven cleanup rc/temp contract failed" >&2
        return 1
    fi
    for state in cleanup-xargs cleanup-watchdog; do
        leader=$(<"$selftest_tmp/$state.pid")
        if process_is_live "$leader"; then
            echo "launcher selftest signal cleanup left $state alive" >&2
            return 1
        fi
    done
    launcher_selftest_cleanup
    trap - EXIT
    printf 'launcher selftest: ok (rcs, cooperative/stubborn groups, signal cleanup)\n'
}

if [[ $# -eq 1 && $1 == --selftest ]]; then
    launcher_selftest
    exit 0
fi

if [[ $# -ne 1 ]]; then
    echo "usage: $0 CAMPAIGN_ROOT" >&2
    exit 2
fi

root=$(realpath -e -- "$1")
launcher=$(realpath -e -- "$0")
prepare=$root/freeze/tools/wh2_h15_v5_holdout_prepare.py
analyzer=$root/freeze/tools/wh2_h15_v5_holdout_analyze.py
runner=$root/freeze/tools/wh2_h15_v5_holdout_job.sh
jobs=$root/recovery/jobs.tsv
workers=120
expected_jobs=17424
expected_artifacts=87123
thermal_csv=${WH2_THERMAL_CSV:-/tmp/wirehair-enoq-thermal.csv}
launch_tmp=
xargs_pid=
xargs_pgid=
watchdog_pid=
cleanup_checks=50
cleanup_delay=0.1
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

launcher_nice=$(ps -o ni= -p $$ | tr -d '[:space:]')
if [[ $launcher_nice != 0 ]]; then
    echo "holdout launcher must run at normal nice level 0 (found $launcher_nice)" >&2
    exit 2
fi
for directory in "$root" "$root/freeze" "$root/recovery" "$root/meta"; do
    if [[ ! -d $directory || -L $directory || $(realpath -e -- "$directory") != "$directory" ]]; then
        echo "campaign path must be a direct nonsymlink directory: $directory" >&2
        exit 2
    fi
done
for file in "$prepare" "$analyzer" "$runner" "$launcher" "$jobs"; do
    if [[ ! -f $file || -L $file ]]; then
        echo "missing or symlinked sealed launcher input: $file" >&2
        exit 2
    fi
done
if [[ $(realpath -e -- "$launcher") != $(realpath -e -- "$root/freeze/tools/wh2_h15_v5_holdout_launch.sh") ]]; then
    echo "only the sealed launcher snapshot may run the holdout" >&2
    exit 2
fi
python3 "$analyzer" audit-inputs --root "$root"
if [[ $(wc -l <"$jobs") -ne $((expected_jobs + 1)) ]]; then
    echo "recovery job ledger count mismatch" >&2
    exit 2
fi

for directory in recovery/raw recovery/status recovery/commands analysis; do
    path=$root/$directory
    if [[ -e $path || -L $path ]]; then
        echo "refusing an existing campaign output directory: $path" >&2
        exit 2
    fi
    mkdir -m 700 -- "$path"
done
for relative in recovery/run_start.json recovery/run_finish.json \
        recovery/thermal_interval.csv recovery/artifact_manifest.sha256 \
        meta/completion_seal.sha256; do
    if [[ -e $root/$relative || -L $root/$relative ]]; then
        echo "refusing to overwrite run artifact: $root/$relative" >&2
        exit 2
    fi
done
if [[ ! -f $thermal_csv || -L $thermal_csv ]]; then
    echo "required live CPU/DIMM thermal CSV is missing or symlinked: $thermal_csv" >&2
    exit 2
fi

launch_tmp=$(mktemp -d -- "$root/.holdout-launch.XXXXXX")
chmod 700 "$launch_tmp"
export PYTHONPYCACHEPREFIX=$launch_tmp/pycache
done_file=$launch_tmp/workers.done
watchdog_output=$launch_tmp/thermal_interval.csv
watchdog_report=$launch_tmp/watchdog.json
thermal_start=$(wc -l <"$thermal_csv")
if ((thermal_start < 2)); then
    echo "thermal CSV does not contain a complete live sample" >&2
    exit 2
fi

python3 "$analyzer" run-record --root "$root" --output "$launch_tmp/run_start.json" \
    --kind start --workers "$workers" --thermal-source "$thermal_csv" \
    --thermal-line "$thermal_start" --launcher-nice "$launcher_nice"
ln -T -- "$launch_tmp/run_start.json" "$root/recovery/run_start.json"
rm -f -- "$launch_tmp/run_start.json"

python3 "$analyzer" watchdog --thermal-source "$thermal_csv" \
    --start-line "$thermal_start" --done-file "$done_file" \
    --output "$watchdog_output" --report "$watchdog_report" &
watchdog_pid=$!

# Positional parameters intentionally expand in the isolated child shell.
# shellcheck disable=SC2016
setsid bash -c '
    set -o pipefail
    tail -n +2 -- "$1" | tr "\n" "\0" | xargs -0 -P 120 -n 1 -- "$2" "$3"
' _ "$jobs" "$runner" "$root" >"$launch_tmp/xargs.stdout" \
    2>"$launch_tmp/xargs.stderr" &
xargs_pid=$!
xargs_pgid=$xargs_pid

watchdog_failed=0
while kill -0 "$xargs_pid" 2>/dev/null; do
    if ! kill -0 "$watchdog_pid" 2>/dev/null; then
        set +e
        wait "$watchdog_pid"
        watchdog_rc=$?
        set -e
        watchdog_pid=
        if ((watchdog_rc != 0)); then
            watchdog_failed=1
            echo "live thermal/EDAC watchdog aborted the recovery campaign" >&2
            terminate_process_group "$xargs_pgid"
            break
        fi
        echo "thermal watchdog ended before recovery workers" >&2
        watchdog_failed=1
        terminate_process_group "$xargs_pgid"
        break
    fi
    sleep 0.5
done
wait_for_leader "$xargs_pid"
xargs_rc=$leader_wait_rc
xargs_pid=
xargs_topology_failed=0
if process_group_alive "$xargs_pgid"; then
    echo "xargs leader exited while recovery process-group members survived" >&2
    terminate_process_group "$xargs_pgid"
    xargs_topology_failed=1
fi
xargs_pgid=
: >"$done_file"
if [[ -n ${watchdog_pid:-} ]]; then
    set +e
    wait "$watchdog_pid"
    watchdog_rc=$?
    set -e
    watchdog_pid=
else
    watchdog_rc=${watchdog_rc:-1}
fi
if ((watchdog_failed != 0 || watchdog_rc != 0)); then
    echo "recovery invalidated by live environment watchdog" >&2
    exit 1
fi
if ((xargs_topology_failed != 0)); then
    echo "recovery invalidated by surviving xargs process-group members" >&2
    exit 1
fi
if [[ -s $launch_tmp/xargs.stdout || -s $launch_tmp/xargs.stderr ]]; then
    echo "xargs/worker orchestration emitted unexpected unsealed output" >&2
    exit 1
fi
ln -T -- "$watchdog_output" "$root/recovery/thermal_interval.csv"
rm -f -- "$watchdog_output"

python3 "$analyzer" run-record --root "$root" --output "$launch_tmp/run_finish.json" \
    --kind finish --xargs-rc "$xargs_rc" --watchdog-report "$watchdog_report"
ln -T -- "$launch_tmp/run_finish.json" "$root/recovery/run_finish.json"
rm -f -- "$launch_tmp/run_finish.json"

python3 "$analyzer" audit-inputs --root "$root"
python3 "$analyzer" audit-artifacts --root "$root"
if ((xargs_rc != 0)); then
    echo "one or more sealed recovery workers failed" >&2
    exit 1
fi
python3 "$analyzer" seal-artifacts --root "$root"
if [[ $(wc -l <"$root/recovery/artifact_manifest.sha256") -ne $((expected_artifacts + 1)) ]]; then
    echo "recovery artifact manifest count mismatch" >&2
    exit 1
fi
python3 "$analyzer" analyze --root "$root"
python3 "$analyzer" audit-inputs --root "$root"
python3 "$analyzer" verify-artifacts --root "$root"
unset PYTHONPYCACHEPREFIX
rm -rf -- "$launch_tmp"
launch_tmp=
python3 "$prepare" seal-completion --root "$root"
printf 'complete\n'
