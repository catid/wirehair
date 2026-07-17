#!/usr/bin/env bash
# Run one arm of one sealed normalized-H15-v5 holdout invocation.
set -euo pipefail
export LC_ALL=C TZ=UTC PYTHONDONTWRITEBYTECODE=1

if [[ $# -ne 2 ]]; then
    echo "usage: $0 CAMPAIGN_ROOT CANONICAL_JOB_ROW" >&2
    exit 2
fi

root=$(realpath -e -- "$1")
line=$2

jobs=$root/recovery/jobs.tsv
binary=$root/freeze/wirehair_v2_bench
analyzer=$root/freeze/tools/wh2_h15_v5_holdout_analyze.py
raw_dir=$root/recovery/raw
status_dir=$root/recovery/status
command_dir=$root/recovery/commands

for directory in "$root" "$root/recovery" "$raw_dir" "$status_dir" "$command_dir"; do
    if [[ ! -d $directory || -L $directory || $(realpath -e -- "$directory") != "$directory" ]]; then
        echo "campaign path must be a direct nonsymlink directory: $directory" >&2
        exit 2
    fi
done
if [[ ! -f $jobs || -L $jobs || ! -x $binary || -L $binary || ! -f $analyzer || -L $analyzer ]]; then
    echo "missing or symlinked sealed job ledger, binary, or analyzer" >&2
    exit 2
fi
IFS=$'\t' read -r manifest_id stem phase group_id seed_index seed schedule loss \
    arm profile trials k_count k_csv extra <<<"$line"
job_id=$manifest_id
if [[ ! $job_id =~ ^(0|[1-9][0-9]{0,4})$ ]] || ((10#$job_id >= 17424)) || \
      [[ -n ${extra:-} || \
      ! $stem =~ ^recovery[0-9]{5}_[A-Za-z0-9_-]+$ || \
      ! $seed =~ ^0x[0-9a-f]{16}$ || ! $seed_index =~ ^[0-2]$ || \
      ! $loss =~ ^0\.(10|35|50|65)$ || \
      ! $schedule =~ ^(iid|burst|permutation|systematic-first|repair-only|adversarial)$ || \
      ! $k_count =~ ^[1-9][0-9]*$ || ! $k_csv =~ ^[0-9]+(,[0-9]+)*$ ]]; then
    echo "malformed recovery job row: $job_id" >&2
    exit 2
fi
if [[ $(awk -F, '{print NF}' <<<"$k_csv") -ne $k_count ]]; then
    echo "recovery job K count mismatch: $job_id" >&2
    exit 2
fi
case $phase:$group_id:$trials in
    breadth:[0-9]*:1)
        if ((10#$group_id >= 120)); then
            echo "invalid breadth group: $job_id" >&2
            exit 2
        fi
        ;;
    depth:changed98:16)
        if [[ $k_count != 98 ]]; then
            echo "invalid depth K count: $job_id" >&2
            exit 2
        fi
        ;;
    *)
        echo "invalid recovery phase/group/trials: $job_id" >&2
        exit 2
        ;;
esac
case $arm:$profile in
    v4:normalized-h15-v4|v5:normalized-h15-v5) ;;
    *)
        echo "invalid arm/profile pair: $job_id" >&2
        exit 2
        ;;
esac

stdout=$raw_dir/$stem.csv
stderr=$raw_dir/$stem.stderr
timing=$raw_dir/$stem.time
status=$status_dir/$stem.ok
command_file=$command_dir/$stem.txt
for path in "$stdout" "$stderr" "$timing" "$status" "$command_file"; do
    if [[ -e $path || -L $path ]]; then
        echo "refusing to overwrite recovery artifact: $path" >&2
        exit 2
    fi
done

tmp_dir=$(mktemp -d -- "$raw_dir/.${stem}.XXXXXX")
chmod 700 "$tmp_dir"
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
    local source=$1 destination=$2
    if ! ln -T -- "$source" "$destination"; then
        echo "concurrent recovery artifact publication refused: $destination" >&2
        exit 2
    fi
    rm -f -- "$source"
}

command=(
    "$binary" precodefail
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
    --packet-peel-seed-table "$profile"
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
    printf 'job_id=%s\tphase=%s\tarm=%s\tseed_index=%s\tk_count=%s\tcommand=' \
        "$job_id" "$phase" "$arm" "$seed_index" "$k_count"
    printf '%q ' "${command[@]}"
    printf '\n'
} >"$tmp_command"

/usr/bin/time -f 'elapsed_s=%e user_s=%U sys_s=%S cpu_pct=%P max_rss_kb=%M exit=%x' \
    -o "$tmp_timing" -- "${command[@]}" >"$tmp_stdout" 2>"$tmp_stderr"
python3 "$analyzer" validate-temp --job-row "$line" \
    --stdout "$tmp_stdout" --stderr "$tmp_stderr"
publish_exclusive "$tmp_stdout" "$stdout"
publish_exclusive "$tmp_stderr" "$stderr"
publish_exclusive "$tmp_timing" "$timing"
publish_exclusive "$tmp_command" "$command_file"
printf 'ok\n' >"$tmp_status"
publish_exclusive "$tmp_status" "$status"
