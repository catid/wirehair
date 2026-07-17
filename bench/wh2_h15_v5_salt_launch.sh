#!/usr/bin/env bash
# Launch the sealed baseline-only exact-witness screen with a 120x1 worker cap.
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "usage: $0 EXPERIMENT BASELINE_BINARY AUDITOR JOB_RUNNER" >&2
    exit 2
fi

experiment=$1
baseline=$2
auditor=$3
job_runner=$4
expected_baseline_sha=4340fd9444b6dc2adbd94304b55e987d9887efd8f6d86f205511e7fda75cb41c
workers=120
temp_artifacts=()
cleanup() {
    local path
    for path in "${temp_artifacts[@]}"; do
        rm -f "$path"
    done
}
trap cleanup EXIT

launcher_nice=$(ps -o ni= -p $$ | tr -d ' ')
if [[ $launcher_nice != 0 ]]; then
    echo "discovery launcher must run at normal nice level 0 (found $launcher_nice)" >&2
    exit 2
fi

for path in "$baseline" "$auditor" "$job_runner"; do
    if [[ ! -f $path ]]; then
        echo "missing required file: $path" >&2
        exit 2
    fi
done
if [[ $(sha256sum "$baseline" | awk '{print $1}') != "$expected_baseline_sha" ]]; then
    echo "baseline binary hash mismatch" >&2
    exit 2
fi

seal=$experiment/meta/discovery_seal.sha256
manifest=$experiment/meta/discovery_manifest.txt
if [[ ! -f $seal || ! -f $manifest ]]; then
    echo "missing discovery seal or manifest" >&2
    exit 2
fi
(cd / && sha256sum -c "$seal") >/dev/null

manifest_value() {
    local key=$1
    awk -F= -v key="$key" '$1 == key {print substr($0, length(key) + 2)}' "$manifest"
}
expected_auditor_path=$(manifest_value auditor_path)
expected_auditor_sha=$(manifest_value auditor_sha256)
expected_job_runner_path=$(manifest_value job_runner_path)
expected_job_runner_sha=$(manifest_value job_runner_sha256)
if [[ $(realpath "$auditor") != "$expected_auditor_path" || \
      $(sha256sum "$auditor" | awk '{print $1}') != "$expected_auditor_sha" || \
      $(realpath "$job_runner") != "$expected_job_runner_path" || \
      $(sha256sum "$job_runner" | awk '{print $1}') != "$expected_job_runner_sha" ]]; then
    echo "supplied auditor/job runner is not the sealed implementation" >&2
    exit 2
fi

jobs=$experiment/meta/exact_jobs.tsv
job_count=$(($(wc -l <"$jobs") - 1))
if [[ $job_count -ne 2304 ]]; then
    echo "expected 2304 jobs, found $job_count" >&2
    exit 2
fi

mkdir -p \
    "$experiment/raw/exact" \
    "$experiment/meta/exact_status" \
    "$experiment/meta/exact_commands" \
    "$experiment/logs"
if find "$experiment/raw/exact" "$experiment/meta/exact_status" \
        "$experiment/meta/exact_commands" -mindepth 1 -print -quit | grep -q .; then
    echo "exact artifact directories are not empty" >&2
    exit 2
fi
run_start=$experiment/meta/run_start.txt
run_finish=$experiment/meta/run_finish.txt
thermal_interval=$experiment/thermal_interval.csv
exit_file=$experiment/exit_code
for path in "$run_start" "$run_finish" "$thermal_interval" "$exit_file"; do
    if [[ -e $path ]]; then
        echo "refusing to overwrite run artifact: $path" >&2
        exit 2
    fi
done

publish_exclusive() {
    local source=$1
    local destination=$2
    if ! ln "$source" "$destination"; then
        echo "concurrent run artifact publication refused: $destination" >&2
        exit 2
    fi
    rm -f "$source"
}

thermal_csv=/tmp/wirehair-enoq-thermal.csv
if [[ ! -f $thermal_csv ]]; then
    echo "required CPU/DIMM thermal worker CSV is missing: $thermal_csv" >&2
    exit 2
fi
thermal_start=$(wc -l <"$thermal_csv")
edac_ce_start=$(awk '{s+=$1} END{print s+0}' /sys/devices/system/edac/mc/mc*/ce_count 2>/dev/null || true)
edac_ue_start=$(awk '{s+=$1} END{print s+0}' /sys/devices/system/edac/mc/mc*/ue_count 2>/dev/null || true)
started=$(date -u +%FT%TZ)

tmp_run_start=$experiment/meta/.run_start.$$.txt
temp_artifacts+=("$tmp_run_start")
{
    printf 'started_utc=%s\n' "$started"
    printf 'workers=%s\nthreads_per_invocation=1\nlauncher_nice=%s\njob_count=%s\n' \
        "$workers" "$launcher_nice" "$job_count"
    printf 'thermal_csv=%s\nthermal_line_start=%s\n' "$thermal_csv" "$thermal_start"
    printf 'edac_ce_start=%s\nedac_ue_start=%s\n' "$edac_ce_start" "$edac_ue_start"
} >"$tmp_run_start"
publish_exclusive "$tmp_run_start" "$run_start"

set +e
seq 0 $((job_count - 1)) |
    xargs -P "$workers" -n 1 -- "$job_runner" "$experiment" "$baseline" "$auditor"
xargs_rc=$?
set -e

thermal_end=$(wc -l <"$thermal_csv")
tmp_thermal=$experiment/.thermal_interval.$$.csv
temp_artifacts+=("$tmp_thermal")
sed -n '1p' "$thermal_csv" >"$tmp_thermal"
sed -n "$((thermal_start + 1)),${thermal_end}p" "$thermal_csv" >>"$tmp_thermal"
publish_exclusive "$tmp_thermal" "$thermal_interval"
edac_ce_end=$(awk '{s+=$1} END{print s+0}' /sys/devices/system/edac/mc/mc*/ce_count 2>/dev/null || true)
edac_ue_end=$(awk '{s+=$1} END{print s+0}' /sys/devices/system/edac/mc/mc*/ue_count 2>/dev/null || true)
finished=$(date -u +%FT%TZ)

status_count=$(find "$experiment/meta/exact_status" -type f -name '*.ok' | wc -l)
stdout_count=$(find "$experiment/raw/exact" -type f -name '*.csv' | wc -l)
stderr_nonempty=$(find "$experiment/raw/exact" -type f -name '*.stderr' -size +0c | wc -l)
stderr_count=$(find "$experiment/raw/exact" -type f -name '*.stderr' | wc -l)
time_count=$(find "$experiment/raw/exact" -type f -name '*.time' | wc -l)
command_count=$(find "$experiment/meta/exact_commands" -type f -name '*.txt' | wc -l)

tmp_run_finish=$experiment/meta/.run_finish.$$.txt
temp_artifacts+=("$tmp_run_finish")
{
    printf 'finished_utc=%s\nxargs_rc=%s\n' "$finished" "$xargs_rc"
    printf 'status_count=%s\nstdout_count=%s\nstderr_count=%s\nstderr_nonempty=%s\n' \
        "$status_count" "$stdout_count" "$stderr_count" "$stderr_nonempty"
    printf 'time_count=%s\ncommand_count=%s\n' "$time_count" "$command_count"
    printf 'thermal_line_end=%s\nedac_ce_end=%s\nedac_ue_end=%s\n' \
        "$thermal_end" "$edac_ce_end" "$edac_ue_end"
} >"$tmp_run_finish"
publish_exclusive "$tmp_run_finish" "$run_finish"

if [[ $xargs_rc -ne 0 || $status_count -ne $job_count || \
      $stdout_count -ne $job_count || $stderr_count -ne $job_count || \
      $time_count -ne $job_count || $command_count -ne $job_count || \
      $stderr_nonempty -ne 0 ]]; then
    echo "discovery worker or post-run artifact count failure" >&2
    exit 1
fi

python3 "$auditor" analyze-exact --experiment "$experiment" --shortlist 16
tmp_exit=$experiment/.exit_code.$$
temp_artifacts+=("$tmp_exit")
printf 'complete\n' >"$tmp_exit"
publish_exclusive "$tmp_exit" "$exit_file"
