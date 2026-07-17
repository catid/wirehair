#!/usr/bin/env bash
# Launch the sealed normalized-H15-v5 selection campaign at 120x1 concurrency.
set -euo pipefail
export LC_ALL=C TZ=UTC PYTHONDONTWRITEBYTECODE=1

if [[ $# -ne 5 ]]; then
    echo "usage: $0 EXPERIMENT BASELINE_BINARY PREPARE ANALYZER JOB_RUNNER" >&2
    exit 2
fi

experiment=$(realpath "$1")
baseline=$(realpath "$2")
prepare=$(realpath "$3")
analyzer=$(realpath "$4")
job_runner=$(realpath "$5")
launcher=$(realpath "$0")
expected_baseline_sha=4340fd9444b6dc2adbd94304b55e987d9887efd8f6d86f205511e7fda75cb41c
expected_header=$'job_id\tstem\tstage\tschedule\tseed_index\tseed\tloss\tarm\tsalt\ttrials\tk_count\tk_csv'
expected_thermal_header='utc,monotonic_s,cpu_busy_pct,cpu_avg_mhz,cpu_tctl_c,dimm_i2c1_50_c,dimm_i2c1_51_c,dimm_i2c1_52_c,dimm_i2c1_53_c,dimm_i2c2_50_c,dimm_i2c2_51_c,dimm_i2c2_52_c,dimm_i2c2_53_c,dimm_read_errors,load1,load5,load15,edac_ce,edac_ue'
workers=120
launch_tmp=
cleanup() {
    if [[ -n $launch_tmp ]]; then
        rm -rf -- "$launch_tmp"
    fi
}
trap cleanup EXIT

launcher_nice=$(ps -o ni= -p $$ | tr -d '[:space:]')
if [[ $launcher_nice != 0 ]]; then
    echo "selection launcher must run at normal nice level 0 (found $launcher_nice)" >&2
    exit 2
fi
if [[ ! -d $experiment || ! -x $baseline || ! -f $prepare || ! -f $analyzer || \
      ! -x $job_runner || ! -f $launcher ]]; then
    echo "missing campaign, baseline, prepare, analyzer, or selection script" >&2
    exit 2
fi
expected_prepare=$(dirname "$analyzer")/wh2_h15_v5_selection_prepare.py
if [[ ! -f $expected_prepare ]] || \
   [[ $(realpath -e -- "$expected_prepare") != "$prepare" ]]; then
    echo "PREPARE must be the analyzer's exact sibling module" >&2
    exit 2
fi
meta=$experiment/meta
if [[ ! -d $meta || -L $meta || $(realpath -e -- "$meta") != "$meta" ]]; then
    echo "campaign meta must be a direct, nonsymlink directory" >&2
    exit 2
fi
if [[ $(sha256sum "$baseline" | awk '{print $1}') != "$expected_baseline_sha" ]]; then
    echo "baseline binary hash mismatch" >&2
    exit 2
fi

seal=$experiment/meta/selection_seal.sha256
if [[ ! -f $seal || -L $seal ]]; then
    echo "missing pre-launch selection seal: $seal" >&2
    exit 2
fi
sealed_files=(
    "$baseline"
    "$prepare"
    "$analyzer"
    "$job_runner"
    "$launcher"
    "$experiment/meta/manifest.json"
    "$experiment/meta/selection_seeds.tsv"
    "$experiment/meta/jobs.tsv"
    "$experiment/meta/groups.tsv"
    "$experiment/meta/discovery_k.csv"
    "$experiment/meta/discovery_failures.csv"
    "$experiment/meta/exact_salt_scores.csv"
    "$experiment/meta/requested_k_salt.tsv"
)
selection_seal_digest=$(sha256sum "$seal" | awk '{print $1}')
verify_selection_seal() {
    local path expected_record
    [[ -f $seal && ! -L $seal && \
       $(sha256sum "$seal" | awk '{print $1}') == "$selection_seal_digest" && \
       $(wc -l <"$seal") -eq ${#sealed_files[@]} ]] || return 1
    for path in "${sealed_files[@]}"; do
        [[ $path == /* && $path != *\\* && $path != *$'\n'* && \
           $path != *$'\r'* && -f $path && ! -L $path ]] || return 1
        expected_record=$(sha256sum "$path") || return 1
        [[ $(grep -Fxc -- "$expected_record" "$seal") -eq 1 ]] || return 1
    done
    (cd / && sha256sum -c "$seal") >/dev/null
}
if ! verify_selection_seal; then
    echo "pre-launch selection seal verification failed" >&2
    exit 2
fi

jobs=$experiment/meta/jobs.tsv
if [[ $(sed -n '1p' "$jobs") != "$expected_header" ]]; then
    echo "jobs.tsv header mismatch" >&2
    exit 2
fi
job_count=$(($(wc -l <"$jobs") - 1))
if [[ $job_count -ne 1980 ]]; then
    echo "expected 1980 jobs, found $job_count" >&2
    exit 2
fi

campaign_dirs=(
    "$experiment/raw"
    "$experiment/meta/status"
    "$experiment/meta/commands"
    "$experiment/logs"
    "$experiment/analysis"
)
for path in "${campaign_dirs[@]}"; do
    if [[ -L $path || -e $path && ! -d $path ]]; then
        echo "campaign output parent must be a direct, nonsymlink directory: $path" >&2
        exit 2
    fi
done
mkdir -p "$experiment/raw" "$experiment/meta/status" \
    "$experiment/meta/commands" "$experiment/logs"
for path in "$experiment/raw" "$experiment/meta/status" \
        "$experiment/meta/commands" "$experiment/logs"; do
    if [[ -L $path || $(realpath -e -- "$path") != "$path" ]]; then
        echo "campaign output parent escaped experiment: $path" >&2
        exit 2
    fi
done
if find "$experiment/raw" "$experiment/meta/status" "$experiment/meta/commands" \
        -mindepth 1 -print -quit | grep -q .; then
    echo "selection artifact directories are not empty" >&2
    exit 2
fi
if [[ -d $experiment/analysis ]] && \
   find "$experiment/analysis" -mindepth 1 -print -quit | grep -q .; then
    echo "campaign analysis directory is not empty" >&2
    exit 2
fi
run_start=$experiment/meta/run_start.txt
run_finish=$experiment/meta/run_finish.txt
thermal_interval=$experiment/thermal_interval.csv
exit_file=$experiment/exit_code
artifact_manifest=$experiment/meta/artifact_manifest.sha256
completion_seal=$experiment/meta/completion_seal.sha256
for path in "$run_start" "$run_finish" "$thermal_interval" "$exit_file" \
        "$artifact_manifest" "$completion_seal"; do
    if [[ -e $path ]]; then
        echo "refusing to overwrite run artifact: $path" >&2
        exit 2
    fi
done
launch_tmp=$(mktemp -d -- "$experiment/.selection-launch.XXXXXX")
export PYTHONPYCACHEPREFIX=$launch_tmp/pycache

publish_exclusive() {
    local source=$1
    local destination=$2
    if ! ln -T -- "$source" "$destination"; then
        echo "concurrent run artifact publication refused: $destination" >&2
        exit 2
    fi
    rm -f "$source"
}

write_sha_manifest() {
    local destination=$1
    shift
    local paths=("$@")
    local offset path
    for path in "${paths[@]}"; do
        [[ $path == /* && $path != *\\* && $path != *$'\n'* && \
           $path != *$'\r'* && -f $path && ! -L $path ]] || return 2
    done
    : >"$destination"
    for ((offset = 0; offset < ${#paths[@]}; offset += 256)); do
        sha256sum "${paths[@]:offset:256}" >>"$destination" || return 1
    done
}

thermal_csv=/tmp/wirehair-enoq-thermal.csv
if [[ ! -f $thermal_csv ]]; then
    echo "required CPU/DIMM thermal worker CSV is missing: $thermal_csv" >&2
    exit 2
fi
thermal_start=$(wc -l <"$thermal_csv")
if ((thermal_start < 1)); then
    echo "thermal CSV is empty: $thermal_csv" >&2
    exit 2
fi
if [[ $(sed -n '1{s/\r$//;p;}' "$thermal_csv") != "$expected_thermal_header" ]]; then
    echo "thermal CSV header/schema mismatch" >&2
    exit 2
fi
shopt -s nullglob
edac_ce_files=(/sys/devices/system/edac/mc/mc*/ce_count)
edac_ue_files=(/sys/devices/system/edac/mc/mc*/ue_count)
shopt -u nullglob
if ((${#edac_ce_files[@]} == 0 || ${#edac_ue_files[@]} == 0)); then
    echo "required EDAC CE/UE counters are unavailable" >&2
    exit 2
fi
edac_sum() {
    local file value total=0
    for file in "$@"; do
        value=$(<"$file") || return 1
        [[ $value =~ ^(0|[1-9][0-9]*)$ ]] || return 1
        total=$((total + 10#$value))
    done
    printf '%s\n' "$total"
}
if ! edac_ce_start=$(edac_sum "${edac_ce_files[@]}") || \
   ! edac_ue_start=$(edac_sum "${edac_ue_files[@]}"); then
    echo "failed to read canonical EDAC CE/UE counters" >&2
    exit 2
fi
started=$(date -u +%FT%TZ)

tmp_run_start=$launch_tmp/run_start.txt
{
    printf 'started_utc=%s\n' "$started"
    printf 'workers=%s\nthreads_per_invocation=1\nlauncher_nice=%s\njob_count=%s\n' \
        "$workers" "$launcher_nice" "$job_count"
    printf 'thermal_csv=%s\nthermal_line_start=%s\n' "$thermal_csv" "$thermal_start"
    printf 'edac_ce_start=%s\nedac_ue_start=%s\n' "$edac_ce_start" "$edac_ue_start"
} >"$tmp_run_start"
publish_exclusive "$tmp_run_start" "$run_start"

run_monotonic_start=$(awk '{print $1}' /proc/uptime)
set +e
seq 0 $((job_count - 1)) |
    xargs -P "$workers" -n 1 -- "$job_runner" "$experiment" "$baseline" "$analyzer"
xargs_rc=$?
set -e

if ! verify_selection_seal; then
    echo "sealed selection inputs changed during worker execution" >&2
    exit 1
fi

run_monotonic_end=$(awk '{print $1}' /proc/uptime)
thermal_end=$(wc -l <"$thermal_csv")
if ((thermal_end <= thermal_start)); then
    echo "thermal sampler did not append a complete sample during selection run" >&2
    exit 1
fi
tmp_thermal=$launch_tmp/thermal_interval.csv
sed -n '1{s/\r$//;p;}' "$thermal_csv" >"$tmp_thermal"
sed -n "$((thermal_start + 1)),${thermal_end}{s/\r$//;p;}" \
    "$thermal_csv" >>"$tmp_thermal"
if ! awk -F, -v start="$run_monotonic_start" -v finish="$run_monotonic_end" '
    function abs(value) { return value < 0 ? -value : value }
    NR == 1 { next }
    NF != 19 || $2 !~ /^[0-9]+([.][0-9]+)?$/ ||
        $3 !~ /^[0-9]+([.][0-9]+)?$/ { exit 1 }
    {
        now = $2 + 0
        busy = $3 + 0
        if (count == 0) {
            first = now
            min_busy = busy
        } else if (now <= last || now - last > 10) {
            exit 2
        }
        if (busy < min_busy) min_busy = busy
        busy_sum += busy
        last = now
        ++count
    }
    END {
        if (count < 2 || abs(first - start) > 10 ||
            abs(finish - last) > 10 || min_busy < 90 ||
            busy_sum / count < 98) exit 3
    }
' "$tmp_thermal"; then
    echo "thermal/load sampling was malformed, discontinuous, or below the full-load gate" >&2
    exit 1
fi
publish_exclusive "$tmp_thermal" "$thermal_interval"
if ! edac_ce_end=$(edac_sum "${edac_ce_files[@]}") || \
   ! edac_ue_end=$(edac_sum "${edac_ue_files[@]}"); then
    echo "failed to reread EDAC CE/UE counters" >&2
    exit 1
fi
finished=$(date -u +%FT%TZ)

status_count=$(find "$experiment/meta/status" -type f -name '*.ok' | wc -l)
stdout_count=$(find "$experiment/raw" -type f -name '*.csv' | wc -l)
stderr_count=$(find "$experiment/raw" -type f -name '*.stderr' | wc -l)
stderr_nonempty=$(find "$experiment/raw" -type f -name '*.stderr' -size +0c | wc -l)
time_count=$(find "$experiment/raw" -type f -name '*.time' | wc -l)
command_count=$(find "$experiment/meta/commands" -type f -name '*.txt' | wc -l)
raw_entry_count=$(find "$experiment/raw" -mindepth 1 -maxdepth 1 -print | wc -l)
status_entry_count=$(find "$experiment/meta/status" -mindepth 1 -maxdepth 1 -print | wc -l)
command_entry_count=$(find "$experiment/meta/commands" -mindepth 1 -maxdepth 1 -print | wc -l)
artifact_names_ok=1
while IFS=$'\t' read -r listed_id stem rest; do
    if [[ ! $listed_id =~ ^[0-9]+$ || ! $stem =~ ^job[0-9]{4}_[A-Za-z0-9_-]+$ || \
          ! -f $experiment/raw/$stem.csv || ! -f $experiment/raw/$stem.stderr || \
          ! -f $experiment/raw/$stem.time || ! -f $experiment/meta/status/$stem.ok || \
          ! -f $experiment/meta/commands/$stem.txt || \
          $(<"$experiment/meta/status/$stem.ok") != ok ]]; then
        artifact_names_ok=0
        break
    fi
done < <(tail -n +2 "$jobs")

tmp_run_finish=$launch_tmp/run_finish.txt
{
    printf 'finished_utc=%s\nxargs_rc=%s\n' "$finished" "$xargs_rc"
    printf 'status_count=%s\nstdout_count=%s\nstderr_count=%s\nstderr_nonempty=%s\n' \
        "$status_count" "$stdout_count" "$stderr_count" "$stderr_nonempty"
    printf 'time_count=%s\ncommand_count=%s\n' "$time_count" "$command_count"
    printf 'raw_entry_count=%s\nstatus_entry_count=%s\ncommand_entry_count=%s\nartifact_names_ok=%s\n' \
        "$raw_entry_count" "$status_entry_count" "$command_entry_count" "$artifact_names_ok"
    printf 'thermal_line_end=%s\nedac_ce_end=%s\nedac_ue_end=%s\n' \
        "$thermal_end" "$edac_ce_end" "$edac_ue_end"
} >"$tmp_run_finish"
publish_exclusive "$tmp_run_finish" "$run_finish"

if [[ $xargs_rc -ne 0 || $status_count -ne $job_count || \
      $stdout_count -ne $job_count || $stderr_count -ne $job_count || \
      $time_count -ne $job_count || $command_count -ne $job_count || \
      $raw_entry_count -ne $((job_count * 3)) || $status_entry_count -ne $job_count || \
      $command_entry_count -ne $job_count || $artifact_names_ok -ne 1 || \
      $stderr_nonempty -ne 0 || $edac_ce_end != "$edac_ce_start" || \
      $edac_ue_end != "$edac_ue_start" ]]; then
    echo "selection worker or post-run artifact count failure" >&2
    exit 1
fi

artifact_files=()
while IFS=$'\t' read -r listed_id stem rest; do
    artifact_files+=(
        "$experiment/raw/$stem.csv"
        "$experiment/raw/$stem.stderr"
        "$experiment/raw/$stem.time"
        "$experiment/meta/status/$stem.ok"
        "$experiment/meta/commands/$stem.txt"
    )
done < <(tail -n +2 "$jobs")
artifact_files+=("$run_start" "$run_finish" "$thermal_interval")
if [[ ${#artifact_files[@]} -ne 9903 ]]; then
    echo "internal artifact manifest count mismatch: ${#artifact_files[@]}" >&2
    exit 1
fi
tmp_artifact_manifest=$launch_tmp/artifact_manifest.sha256
if ! write_sha_manifest "$tmp_artifact_manifest" "${artifact_files[@]}" || \
   [[ $(wc -l <"$tmp_artifact_manifest") -ne 9903 ]] || \
   ! (cd / && sha256sum -c "$tmp_artifact_manifest") >/dev/null; then
    echo "failed to create or verify exact artifact manifest" >&2
    exit 1
fi
publish_exclusive "$tmp_artifact_manifest" "$artifact_manifest"
artifact_manifest_digest=$(sha256sum "$artifact_manifest" | awk '{print $1}')

python3 "$analyzer" analyze --experiment "$experiment"
if [[ $(sha256sum "$artifact_manifest" | awk '{print $1}') != "$artifact_manifest_digest" ]] || \
   ! (cd / && sha256sum -c "$artifact_manifest") >/dev/null; then
    echo "sealed run inputs changed during analysis" >&2
    exit 1
fi
if ! verify_selection_seal; then
    echo "sealed selection inputs changed during analysis" >&2
    exit 1
fi
analysis_dir=$experiment/analysis
analysis_files=(
    "$analysis_dir/decision_ledger.csv"
    "$analysis_dir/selection_table.tsv"
    "$analysis_dir/work_metrics.csv"
    "$analysis_dir/recovery_metrics.csv"
    "$analysis_dir/summary.json"
)
if [[ $(find "$analysis_dir" -mindepth 1 -maxdepth 1 -print | wc -l) -ne 5 ]]; then
    echo "analysis output set must contain exactly five artifacts" >&2
    exit 1
fi

completion_files=("$seal" "$artifact_manifest" "${analysis_files[@]}")
tmp_completion=$launch_tmp/completion_seal.sha256
if [[ ${#completion_files[@]} -ne 7 ]] || \
   ! write_sha_manifest "$tmp_completion" "${completion_files[@]}" || \
   [[ $(wc -l <"$tmp_completion") -ne 7 ]] || \
   ! (cd / && sha256sum -c "$tmp_completion") >/dev/null; then
    echo "failed to create or verify completion seal" >&2
    exit 1
fi
publish_exclusive "$tmp_completion" "$completion_seal"

tmp_exit=$launch_tmp/exit_code
printf 'complete\n' >"$tmp_exit"
publish_exclusive "$tmp_exit" "$exit_file"
