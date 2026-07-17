#!/usr/bin/env bash
# Run one immutable, baseline-only normalized-H15 exact-witness salt probe.
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "usage: $0 EXPERIMENT BASELINE_BINARY AUDITOR JOB_ID" >&2
    exit 2
fi

experiment=$1
baseline=$2
auditor=$3
job_id=$4

if [[ ! $job_id =~ ^[0-9]+$ ]]; then
    echo "invalid job ID: $job_id" >&2
    exit 2
fi

jobs=$experiment/meta/exact_jobs.tsv
if [[ ! -f $jobs || ! -x $baseline || ! -f $auditor ]]; then
    echo "missing jobs, executable baseline, or auditor" >&2
    exit 2
fi

line=$(sed -n "$((job_id + 2))p" "$jobs")
if [[ -z $line ]]; then
    echo "job ID outside manifest: $job_id" >&2
    exit 2
fi
IFS=$'\t' read -r manifest_id stem schedule seed_index seed loss salt k_count k_csv <<<"$line"
if [[ $manifest_id != "$job_id" || -z $stem || -z $k_csv ]]; then
    echo "malformed or mismatched job row: $job_id" >&2
    exit 2
fi

raw_dir=$experiment/raw/exact
status_dir=$experiment/meta/exact_status
command_dir=$experiment/meta/exact_commands
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

tmp_prefix=$raw_dir/.$stem.$$
tmp_stdout=$tmp_prefix.csv
tmp_stderr=$tmp_prefix.stderr
tmp_timing=$tmp_prefix.time
tmp_command=$command_dir/.$stem.$$.txt
tmp_status=$status.$$
cleanup() {
    rm -f "$tmp_stdout" "$tmp_stderr" "$tmp_timing" "$tmp_command" "$tmp_status"
}
trap cleanup EXIT

publish_exclusive() {
    local source=$1
    local destination=$2
    if ! ln "$source" "$destination"; then
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
    --trials 1
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
    printf 'arm=baseline-only\tjob_id=%s\tseed_index=%s\tk_count=%s\tcommand=' \
        "$job_id" "$seed_index" "$k_count"
    printf '%q ' "${command[@]}"
    printf '\n'
} >"$tmp_command"

/usr/bin/time -f 'elapsed_s=%e user_s=%U sys_s=%S cpu_pct=%P max_rss_kb=%M exit=%x' \
    -o "$tmp_timing" -- "${command[@]}" >"$tmp_stdout" 2>"$tmp_stderr"

publish_exclusive "$tmp_stdout" "$stdout"
publish_exclusive "$tmp_stderr" "$stderr"
publish_exclusive "$tmp_timing" "$timing"
publish_exclusive "$tmp_command" "$command_file"

python3 "$auditor" validate-exact-job \
    --experiment "$experiment" --job-id "$job_id"
printf 'ok\n' >"$tmp_status"
publish_exclusive "$tmp_status" "$status"
