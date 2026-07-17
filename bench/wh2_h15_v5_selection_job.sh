#!/usr/bin/env bash
# Run one sealed normalized-H15-v5 Stage G/Stage T selection invocation.
set -euo pipefail
export LC_ALL=C TZ=UTC PYTHONDONTWRITEBYTECODE=1

if [[ $# -ne 4 ]]; then
    echo "usage: $0 EXPERIMENT BASELINE_BINARY ANALYZER JOB_ID" >&2
    exit 2
fi

experiment=$1
baseline=$2
analyzer=$3
job_id=$4
expected_baseline_sha=4340fd9444b6dc2adbd94304b55e987d9887efd8f6d86f205511e7fda75cb41c
expected_header=$'job_id\tstem\tstage\tschedule\tseed_index\tseed\tloss\tarm\tsalt\ttrials\tk_count\tk_csv'

if [[ ! $job_id =~ ^(0|[1-9][0-9]{0,3})$ ]]; then
    echo "invalid job ID: $job_id" >&2
    exit 2
fi
job_number=$job_id
if ((job_number >= 1980)); then
    echo "job ID outside [0,1980): $job_id" >&2
    exit 2
fi

jobs=$experiment/meta/jobs.tsv
seeds=$experiment/meta/selection_seeds.tsv
raw_dir=$experiment/raw
status_dir=$experiment/meta/status
command_dir=$experiment/meta/commands
if [[ ! -f $jobs || ! -f $seeds || ! -x $baseline || ! -f $analyzer || \
      ! -d $raw_dir || ! -d $status_dir || ! -d $command_dir ]]; then
    echo "missing campaign metadata, artifact directories, baseline, or analyzer" >&2
    exit 2
fi
if [[ $(sha256sum "$baseline" | awk '{print $1}') != "$expected_baseline_sha" ]]; then
    echo "baseline binary hash mismatch" >&2
    exit 2
fi
if [[ $(sed -n '1p' "$jobs") != "$expected_header" || $(wc -l <"$jobs") -ne 1981 ]]; then
    echo "jobs.tsv header or invocation count mismatch" >&2
    exit 2
fi

line=$(sed -n "$((job_number + 2))p" "$jobs")
IFS=$'\t' read -r manifest_id stem stage schedule seed_index seed loss arm salt \
    trials k_count k_csv extra <<<"$line"
if [[ $manifest_id != "$job_id" || -n ${extra:-} || -z $stem || -z $k_csv || \
      ! $stem =~ ^job[0-9]{4}_[A-Za-z0-9_-]+$ || \
      ! $seed =~ ^0x[0-9a-f]{16}$ || ! $k_count =~ ^[1-9][0-9]*$ || \
      ! $k_csv =~ ^[0-9]+(,[0-9]+)*$ ]]; then
    echo "malformed or mismatched job row: $job_id" >&2
    exit 2
fi
if [[ $(awk -F, '{print NF}' <<<"$k_csv") -ne $k_count ]]; then
    echo "job K count mismatch: $job_id" >&2
    exit 2
fi

case $stage in
    G)
        if [[ $seed_index != 0 || $loss != 0.50 || $trials != 1 || \
              ! $schedule =~ ^(burst|adversarial|repair-only)$ ]]; then
            echo "invalid Stage G job geometry: $job_id" >&2
            exit 2
        fi
        ;;
    T)
        if [[ ! $seed_index =~ ^[12]$ || ! $loss =~ ^0\.(35|50|65)$ || \
              $trials != 2 || \
              ! $schedule =~ ^(iid|burst|permutation|systematic-first|repair-only|adversarial)$ ]]; then
            echo "invalid Stage T job geometry: $job_id" >&2
            exit 2
        fi
        ;;
    *)
        echo "invalid job stage: $job_id" >&2
        exit 2
        ;;
esac

case $arm:$salt in
    baseline:0|c27:39|c79:121|c6f:111|ca8:168) ;;
    *)
        echo "invalid arm/salt pair: $job_id" >&2
        exit 2
        ;;
esac
seed_from_meta=$(awk -F '\t' -v wanted="$seed_index" \
    '$1 == "selection" && $2 == wanted {if (found++) exit 3; value=$5} END {if (found != 1) exit 4; print value}' \
    "$seeds") || {
    echo "selection seed metadata mismatch: $job_id" >&2
    exit 2
}
if [[ $seed != "$seed_from_meta" ]]; then
    echo "job seed differs from selection metadata: $job_id" >&2
    exit 2
fi

stdout=$raw_dir/$stem.csv
stderr=$raw_dir/$stem.stderr
timing=$raw_dir/$stem.time
status=$status_dir/$stem.ok
command_file=$command_dir/$stem.txt
for path in "$stdout" "$stderr" "$timing" "$status" "$command_file"; do
    if [[ -e $path ]]; then
        echo "refusing to overwrite artifact: $path" >&2
        exit 2
    fi
done

tmp_dir=$(mktemp -d -- "$raw_dir/.${stem}.XXXXXX")
export PYTHONPYCACHEPREFIX=$tmp_dir/pycache
tmp_stdout=$tmp_dir/stdout.csv
tmp_stderr=$tmp_dir/stderr
tmp_timing=$tmp_dir/time
tmp_command=$tmp_dir/command.txt
tmp_status=$tmp_dir/status.ok
cleanup() {
    rm -rf -- "$tmp_dir"
}
trap cleanup EXIT

publish_exclusive() {
    local source=$1
    local destination=$2
    if ! ln -T -- "$source" "$destination"; then
        echo "concurrent artifact publication refused: $destination" >&2
        exit 2
    fi
    rm -f "$source"
}

command=(
    "$baseline" precodefail
    --N "$k_csv"
    --bb-list 64
    --seed-block-bytes 1280
    --overhead 0
    --trials "$trials"
    --threads 1
    --loss "$loss"
    --seed "$seed"
    --schedule "$schedule"
    --completion mixed
    --mix-count 2
    --packet-peel-seed-xor "$salt"
    --mixed-gf256-rows 11
    --mixed-gf16-rows 4
    --mixed-period 32
    --mixed-geometry shared-x
    --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 68
    --mixed-residue-hash-keyed
    --mixed-independent-extension-residues
    --mixed-extension-residue-seed-xor 78
)

{
    printf 'job_id=%s\tstage=%s\tarm=%s\tseed_index=%s\tk_count=%s\tcommand=' \
        "$job_id" "$stage" "$arm" "$seed_index" "$k_count"
    printf '%q ' "${command[@]}"
    printf '\n'
} >"$tmp_command"

/usr/bin/time -f 'elapsed_s=%e user_s=%U sys_s=%S cpu_pct=%P max_rss_kb=%M exit=%x' \
    -o "$tmp_timing" -- "${command[@]}" >"$tmp_stdout" 2>"$tmp_stderr"
python3 "$analyzer" validate-output --experiment "$experiment" --job-id "$job_id" \
    --stdout "$tmp_stdout" --stderr "$tmp_stderr"
publish_exclusive "$tmp_stdout" "$stdout"
publish_exclusive "$tmp_stderr" "$stderr"
publish_exclusive "$tmp_timing" "$timing"
publish_exclusive "$tmp_command" "$command_file"

python3 "$analyzer" validate-job --experiment "$experiment" --job-id "$job_id"
printf 'ok\n' >"$tmp_status"
publish_exclusive "$tmp_status" "$status"
