#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DRIVER="$ROOT_DIR/tables/run_full_regeneration.sh"
PYTHON=${PYTHON:-python3}

CPUSET=${CPUSET:-0-111}
JOBS=${JOBS:-112}
OMP_THREADS=${OMP_THREADS:-1}
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build/jwg-fixup-free"}

DCOUNT_SEEDS=${DCOUNT_SEEDS:-123456789,987654321,424242424,2718281828}
DCOUNT_TRIALS=${DCOUNT_TRIALS:-5000}
DCOUNT_MAX_FAILURES=${DCOUNT_MAX_FAILURES:-100}
DCOUNT_LOW_COUNT_RUN=${DCOUNT_LOW_COUNT_RUN:-4}
SMALL_DENSE_SEED=${SMALL_DENSE_SEED:-3141592653}
SMALL_DENSE_TRIALS=${SMALL_DENSE_TRIALS:-3500}
MOST_DENSE_SEED=${MOST_DENSE_SEED:-1618033988}
MOST_DENSE_TRIALS=${MOST_DENSE_TRIALS:-5000}
PEEL_SEED=${PEEL_SEED:-1414213562}
PEEL_TRIALS=${PEEL_TRIALS:-30}
PEEL_MAX_TRIES=${PEEL_MAX_TRIES:-32}
HEAVY_TRIALS=${HEAVY_TRIALS:-1000000}

export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

die() {
    echo "run_full_regeneration: $*" >&2
    exit 2
}

require_uint() {
    local name=$1 value=$2
    [[ "$value" =~ ^[0-9]+$ ]] || die "$name must be an unsigned integer"
}

require_positive() {
    local name=$1 value=$2
    require_uint "$name" "$value"
    ((value > 0)) || die "$name must be positive"
}

validate_environment() {
    require_positive JOBS "$JOBS"
    require_positive OMP_THREADS "$OMP_THREADS"
    require_positive DCOUNT_TRIALS "$DCOUNT_TRIALS"
    require_uint DCOUNT_MAX_FAILURES "$DCOUNT_MAX_FAILURES"
    require_positive DCOUNT_LOW_COUNT_RUN "$DCOUNT_LOW_COUNT_RUN"
    require_uint SMALL_DENSE_SEED "$SMALL_DENSE_SEED"
    require_positive SMALL_DENSE_TRIALS "$SMALL_DENSE_TRIALS"
    require_uint MOST_DENSE_SEED "$MOST_DENSE_SEED"
    require_positive MOST_DENSE_TRIALS "$MOST_DENSE_TRIALS"
    require_uint PEEL_SEED "$PEEL_SEED"
    require_positive PEEL_TRIALS "$PEEL_TRIALS"
    require_positive PEEL_MAX_TRIES "$PEEL_MAX_TRIES"
    require_uint HEAVY_TRIALS "$HEAVY_TRIALS"
    [[ "$DCOUNT_SEEDS" =~ ^[0-9]+(,[0-9]+){3}$ ]] ||
        die "DCOUNT_SEEDS must contain exactly four decimal seeds"
    [[ "$(tr ',' '\n' <<< "$DCOUNT_SEEDS" | sort -u | wc -l)" -eq 4 ]] ||
        die "DCOUNT_SEEDS must be distinct"
    ((DCOUNT_MAX_FAILURES <= DCOUNT_TRIALS)) ||
        die "DCOUNT_MAX_FAILURES cannot exceed DCOUNT_TRIALS"
    ((PEEL_MAX_TRIES <= 256)) || die "PEEL_MAX_TRIES cannot exceed 256"
    local available
    available=$(taskset -c "$CPUSET" nproc) || die "invalid CPUSET=$CPUSET"
    ((available <= 112)) || die "CPUSET exposes $available CPUs; maximum is 112"
    ((JOBS * OMP_THREADS <= available)) ||
        die "JOBS*OMP_THREADS exceeds CPUSET capacity $available"
    [[ -f "$BUILD_DIR/CMakeCache.txt" ]] || die "missing CMake build: $BUILD_DIR"
    grep -q '^WIREHAIR_DISABLE_SEED_FIXUPS:BOOL=ON$' \
        "$BUILD_DIR/CMakeCache.txt" ||
        die "BUILD_DIR must set WIREHAIR_DISABLE_SEED_FIXUPS=ON"
}

binary_for_phase() {
    case "$1" in
        dense-count) echo "$BUILD_DIR/gen_dcounts" ;;
        small) echo "$BUILD_DIR/gen_small_dseeds" ;;
        dense-seeds) echo "$BUILD_DIR/gen_most_dseeds" ;;
        peel) echo "$BUILD_DIR/gen_peel_seeds" ;;
        heavy) echo "$BUILD_DIR/gen_tables" ;;
        *) die "unknown phase: $1" ;;
    esac
}

kind_for_phase() {
    case "$1" in
        dense-count|small|dense-seeds|peel) echo "$1" ;;
        heavy) echo heavy ;;
        *) die "unknown phase: $1" ;;
    esac
}

result_name_for_phase() {
    case "$1" in
        dense-count) echo dense_count.results.tsv ;;
        small) echo small_dense_seeds.results.tsv ;;
        dense-seeds) echo most_dense_seeds.results.tsv ;;
        peel) echo peel_seeds.results.tsv ;;
        heavy) echo heavy.out ;;
        *) die "unknown phase: $1" ;;
    esac
}

phase_parameters() {
    local phase=$1
    case "$phase" in
        dense-count)
            printf 'seeds=%s\ntrials=%s\nmax_failures=%s\nlow_count_run=%s\n' \
                "$DCOUNT_SEEDS" "$DCOUNT_TRIALS" \
                "$DCOUNT_MAX_FAILURES" "$DCOUNT_LOW_COUNT_RUN"
            ;;
        small)
            printf 'seed=%s\ntrials=%s\nN=2..2047\n' \
                "$SMALL_DENSE_SEED" "$SMALL_DENSE_TRIALS"
            ;;
        dense-seeds)
            printf 'seed=%s\ntrials=%s\nindices=2..99\n' \
                "$MOST_DENSE_SEED" "$MOST_DENSE_TRIALS"
            ;;
        peel)
            printf 'seed=%s\ntrials=%s\nmax_tries=%s\nsubdivisions=0..2047\nN=2048..64000\n' \
                "$PEEL_SEED" "$PEEL_TRIALS" "$PEEL_MAX_TRIES"
            ;;
        heavy)
            printf 'trials=%s\nrepeats=2\n' "$HEAVY_TRIALS"
            ;;
    esac
}

write_phase_manifest() {
    local phase=$1 phase_dir=$2 binary=$3
    local expected="$phase_dir/.manifest.expected.$$"
    {
        echo '# wirehair full table regeneration phase v1'
        echo "phase=$phase"
        echo "git_commit=$(git -C "$ROOT_DIR" rev-parse HEAD)"
        echo "cpuset=$CPUSET"
        echo "jobs=$JOBS"
        echo "omp_threads=$OMP_THREADS"
        echo "build_dir=$BUILD_DIR"
        echo "binary=$binary"
        echo "binary_sha256=$(sha256sum "$binary" | awk '{print $1}')"
        echo "cmake_cache_sha256=$(sha256sum "$BUILD_DIR/CMakeCache.txt" | awk '{print $1}')"
        echo "driver_sha256=$(sha256sum "$DRIVER" | awk '{print $1}')"
        echo "validator_sha256=$(sha256sum "$ROOT_DIR/tables/regeneration_results.py" | awk '{print $1}')"
        echo "wirehair_tools_sha256=$(sha256sum "$ROOT_DIR/WirehairTools.cpp" | awk '{print $1}')"
        echo "wirehair_codec_sha256=$(sha256sum "$ROOT_DIR/WirehairCodec.cpp" | awk '{print $1}')"
        case "$phase" in
            dense-count) sha256sum "$ROOT_DIR/tables/GenerateDenseCount.cpp" ;;
            small) sha256sum "$ROOT_DIR/tables/GenerateSmallDenseSeeds.cpp" ;;
            dense-seeds) sha256sum "$ROOT_DIR/tables/GenerateMostDenseSeeds.cpp" ;;
            peel) sha256sum "$ROOT_DIR/tables/GeneratePeelSeeds.cpp" ;;
            heavy) sha256sum "$ROOT_DIR/tables/TableGenerator.cpp" "$ROOT_DIR/tables/HeavyRowGenerator.cpp" ;;
        esac | sed 's/^/source_sha256=/'
        phase_parameters "$phase"
    } > "$expected"
    if [[ -f "$phase_dir/manifest.txt" ]]; then
        cmp -s "$expected" "$phase_dir/manifest.txt" || {
            diff -u "$phase_dir/manifest.txt" "$expected" >&2 || true
            rm -f "$expected"
            die "phase provenance changed; use a new campaign directory"
        }
        rm -f "$expected"
    else
        mv -f "$expected" "$phase_dir/manifest.txt"
    fi
}

verify_complete_shard() {
    local directory=$1
    [[ -f "$directory/COMPLETE" && -f "$directory/checksums.sha256" ]] ||
        return 1
    (cd "$directory" && sha256sum -c checksums.sha256 >/dev/null)
}

run_worker() {
    local phase=$1 key=$2 campaign_dir=$3
    validate_environment
    local binary phase_dir shard_root complete attempt result_name kind
    binary=$(binary_for_phase "$phase")
    phase_dir="$campaign_dir/$phase"
    shard_root="$phase_dir/shards"
    complete="$shard_root/$key.complete"
    result_name=$(result_name_for_phase "$phase")
    kind=$(kind_for_phase "$phase")
    if [[ -d "$complete" ]]; then
        verify_complete_shard "$complete" || die "corrupt completed shard: $complete"
        return 0
    fi
    attempt="$shard_root/.$key.attempt.$BASHPID"
    mkdir -p "$attempt"
    local cpus=$CPUSET threads=$OMP_THREADS
    local -a command
    case "$phase" in
        dense-count)
            local seed=${key#seed-}
            seed=${seed#repeat-}
            IFS=',' read -r -a seed_list <<< "$DCOUNT_SEEDS"
            local position=-1
            for i in "${!seed_list[@]}"; do
                if [[ "${seed_list[i]}" == "$seed" ]]; then position=$i; break; fi
            done
            ((position >= 0)) || die "unknown dense-count seed $seed"
            [[ "$CPUSET" =~ ^([0-9]+)-([0-9]+)$ ]] ||
                die "dense-count requires a contiguous CPUSET range"
            local cpu_first=${BASH_REMATCH[1]} cpu_last=${BASH_REMATCH[2]}
            local cpu_count=$((cpu_last - cpu_first + 1))
            ((cpu_count % 4 == 0)) ||
                die "dense-count CPUSET size must be divisible by four"
            local width=$((cpu_count / 4))
            local lo=$((cpu_first + position * width))
            local hi=$((lo + width - 1))
            cpus="$lo-$hi"
            threads=$width
            if [[ "$key" == repeat-* ]]; then
                cpus=$CPUSET
                threads=$cpu_count
            fi
            command=("$binary" --seed "$seed" --trials "$DCOUNT_TRIALS"
                --nlo 2 --nhi 64000 --max-failures "$DCOUNT_MAX_FAILURES"
                --low-count-run "$DCOUNT_LOW_COUNT_RUN"
                --results "$attempt/$result_name")
            ;;
        small)
            command=("$binary" --seed "$SMALL_DENSE_SEED"
                --trials "$SMALL_DENSE_TRIALS" --nlo "$key" --nhi "$key"
                --nlist "$key" --results "$attempt/$result_name")
            ;;
        dense-seeds)
            command=("$binary" --seed "$MOST_DENSE_SEED"
                --trials "$MOST_DENSE_TRIALS" --dense-index-lo "$key"
                --dense-index-hi "$key" --results "$attempt/$result_name")
            ;;
        peel)
            command=("$binary" --seed "$PEEL_SEED" --trials "$PEEL_TRIALS"
                --sublo "$key" --subhi "$key" --nlo 2048 --nhi 64000
                --max-tries "$PEEL_MAX_TRIES"
                --results "$attempt/$result_name")
            ;;
        heavy)
            command=("$binary" --no-benchmarks --heavy-trials "$HEAVY_TRIALS")
            ;;
    esac
    {
        echo "utc_start=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "cpuset=$cpus"
        echo "omp_threads=$threads"
        printf 'command='
        printf '%q ' "${command[@]}"
        printf '\n'
    } > "$attempt/run.txt"
    if ! taskset -c "$cpus" env OMP_NUM_THREADS="$threads" \
        "${command[@]}" > "$attempt/$phase.out" 2> "$attempt/$phase.err"
    then
        echo "exit=failure" >> "$attempt/run.txt"
        return 1
    fi
    if [[ "$phase" != heavy ]]; then
        "$PYTHON" "$ROOT_DIR/tables/regeneration_results.py" "$kind" \
            "$attempt/$result_name" > /dev/null
    fi
    echo "utc_end=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$attempt/run.txt"
    (
        cd "$attempt"
        files=(run.txt "$phase.out" "$phase.err")
        if [[ "$result_name" != "$phase.out" ]]; then
            files+=("$result_name")
        fi
        sha256sum "${files[@]}" > checksums.sha256
        printf 'complete\n' > COMPLETE
    )
    if [[ -e "$complete" ]]; then
        die "completed shard appeared concurrently: $complete"
    fi
    mv -f "$attempt" "$complete"
    verify_complete_shard "$complete" || die "failed post-promotion verification"
}

canonicalize_phase() {
    local phase=$1 campaign_dir=$2 expected_count=$3
    local phase_dir="$campaign_dir/$phase" kind result_name
    kind=$(kind_for_phase "$phase")
    result_name=$(result_name_for_phase "$phase")
    mapfile -t directories < <(
        if [[ "$phase" == dense-count ]]; then
            find "$phase_dir/shards" -mindepth 1 -maxdepth 1 \
                -type d -name 'seed-*.complete' -print
        else
            find "$phase_dir/shards" -mindepth 1 -maxdepth 1 \
                -type d -name '*.complete' -print
        fi | LC_ALL=C sort
    )
    [[ "${#directories[@]}" -eq "$expected_count" ]] ||
        die "$phase completed ${#directories[@]} shards; expected $expected_count"
    local directory
    for directory in "${directories[@]}"; do
        verify_complete_shard "$directory" || die "corrupt shard: $directory"
    done
    local temporary="$phase_dir/.merged.results.tsv.$$"
    local -a options=(--full --output "$temporary")
    if [[ "$phase" == dense-count ]]; then
        options+=(--min-runs 4)
    fi
    "$PYTHON" "$ROOT_DIR/tables/regeneration_results.py" "$kind" \
        "${directories[@]}" "${options[@]}"
    if [[ -f "$phase_dir/merged.results.tsv" ]]; then
        cmp -s "$temporary" "$phase_dir/merged.results.tsv" ||
            die "$phase canonical output changed"
        rm -f "$temporary"
    else
        mv -f "$temporary" "$phase_dir/merged.results.tsv"
    fi
    sha256sum "$phase_dir/merged.results.tsv" > "$phase_dir/merged.sha256"
    printf 'complete %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
        > "$phase_dir/COMPLETE.tmp"
    mv -f "$phase_dir/COMPLETE.tmp" "$phase_dir/COMPLETE"
}

fingerprint_value() {
    local name=$1
    "$PYTHON" "$ROOT_DIR/tables/apply_regeneration_results.py" fingerprint \
        --source "$ROOT_DIR/WirehairTools.cpp" |
        sed -n "s/^${name}=//p"
}

check_phase_dependency() {
    local phase=$1 campaign_dir=$2 marker name expected actual
    case "$phase" in
        dense-seeds)
            marker="$campaign_dir/dense-count/APPLIED"
            name=dense_graph_hash
            ;;
        peel)
            marker="$campaign_dir/dense-seeds/APPLIED"
            name=dense_seed_hash
            ;;
        *) return 0 ;;
    esac
    [[ -f "$marker" ]] ||
        die "$phase requires the prior training table to be applied: $marker"
    expected=$(sed -n "s/^${name}=//p" "$marker")
    actual=$(fingerprint_value "$name")
    [[ -n "$expected" && "$actual" == "$expected" ]] ||
        die "$phase source fingerprint does not match prior APPLIED marker"
}

mark_applied() {
    local phase=$1 campaign_dir=$2 phase_dir marker expected
    validate_environment
    campaign_dir=$(cd "$campaign_dir" && pwd)
    phase_dir="$campaign_dir/$phase"
    [[ -f "$phase_dir/COMPLETE" ]] || die "$phase is not complete"
    case "$phase" in
        dense-count)
            cmake --build "$BUILD_DIR" --target gen_most_dseeds -j 16
            ;;
        dense-seeds)
            cmake --build "$BUILD_DIR" --target gen_peel_seeds -j 16
            ;;
        small|peel) ;;
        *) die "cannot mark phase applied: $phase" ;;
    esac
    marker="$phase_dir/APPLIED"
    expected="$phase_dir/.APPLIED.expected.$$"
    {
        echo '# provisional training-table application; final holdout still required'
        echo "phase=$phase"
        echo "utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "merged_sha256=$(sha256sum "$phase_dir/merged.results.tsv" 2>/dev/null | awk '{print $1}')"
        echo "wirehair_tools_sha256=$(sha256sum "$ROOT_DIR/WirehairTools.cpp" | awk '{print $1}')"
        "$PYTHON" "$ROOT_DIR/tables/apply_regeneration_results.py" fingerprint \
            --source "$ROOT_DIR/WirehairTools.cpp"
        if [[ "$phase" == dense-count ]]; then
            echo "consumer_binary_sha256=$(sha256sum "$BUILD_DIR/gen_most_dseeds" | awk '{print $1}')"
        elif [[ "$phase" == dense-seeds ]]; then
            echo "consumer_binary_sha256=$(sha256sum "$BUILD_DIR/gen_peel_seeds" | awk '{print $1}')"
        fi
    } > "$expected"
    if [[ -f "$marker" ]]; then
        # Ignore the timestamp when confirming an existing marker.
        diff -u <(grep -v '^utc=' "$marker") \
            <(grep -v '^utc=' "$expected") >/dev/null ||
            die "$phase APPLIED marker does not match current tables"
        rm -f "$expected"
    else
        mv -f "$expected" "$marker"
    fi
    echo "$phase provisional application recorded: $marker"
}

run_phase() {
    local phase=$1 campaign_dir=$2 binary phase_dir
    validate_environment
    binary=$(binary_for_phase "$phase")
    [[ -x "$binary" ]] || die "missing executable: $binary"
    mkdir -p "$campaign_dir" "$campaign_dir/$phase/shards"
    campaign_dir=$(cd "$campaign_dir" && pwd)
    phase_dir="$campaign_dir/$phase"
    check_phase_dependency "$phase" "$campaign_dir"
    write_phase_manifest "$phase" "$phase_dir" "$binary"
    if [[ -f "$phase_dir/COMPLETE" ]]; then
        echo "$phase already complete: $phase_dir"
        return 0
    fi
    export ROOT_DIR DRIVER PYTHON CPUSET JOBS OMP_THREADS BUILD_DIR
    export DCOUNT_SEEDS DCOUNT_TRIALS DCOUNT_MAX_FAILURES DCOUNT_LOW_COUNT_RUN
    export SMALL_DENSE_SEED SMALL_DENSE_TRIALS MOST_DENSE_SEED MOST_DENSE_TRIALS
    export PEEL_SEED PEEL_TRIALS PEEL_MAX_TRIES HEAVY_TRIALS
    case "$phase" in
        dense-count)
            IFS=',' read -r -a seeds <<< "$DCOUNT_SEEDS"
            printf 'seed-%s\n' "${seeds[@]}" | xargs -r -n 1 -P 4 \
                "$DRIVER" __worker dense-count "$campaign_dir"
            "$DRIVER" __worker dense-count "$campaign_dir" \
                "repeat-${seeds[0]}"
            cmp -s \
                "$phase_dir/shards/seed-${seeds[0]}.complete/dense_count.results.tsv" \
                "$phase_dir/shards/repeat-${seeds[0]}.complete/dense_count.results.tsv" ||
                die "dense-count deterministic repeat differs"
            canonicalize_phase dense-count "$campaign_dir" 4
            ;;
        small)
            seq 2 2047 | xargs -r -n 1 -P "$JOBS" \
                "$DRIVER" __worker small "$campaign_dir"
            canonicalize_phase small "$campaign_dir" 2046
            ;;
        dense-seeds)
            seq 2 99 | xargs -r -n 1 -P "$JOBS" \
                "$DRIVER" __worker dense-seeds "$campaign_dir"
            canonicalize_phase dense-seeds "$campaign_dir" 98
            ;;
        peel)
            seq 0 2047 | xargs -r -n 1 -P "$JOBS" \
                "$DRIVER" __worker peel "$campaign_dir"
            canonicalize_phase peel "$campaign_dir" 2048
            ;;
        heavy)
            printf 'first\nsecond\n' | xargs -r -n 1 -P 2 \
                "$DRIVER" __worker heavy "$campaign_dir"
            local first="$phase_dir/shards/first.complete/heavy.out"
            local second="$phase_dir/shards/second.complete/heavy.out"
            cmp -s "$first" "$second" || die "heavy generator repeats differ"
            sha256sum "$first" > "$phase_dir/merged.sha256"
            printf 'complete %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
                > "$phase_dir/COMPLETE.tmp"
            mv -f "$phase_dir/COMPLETE.tmp" "$phase_dir/COMPLETE"
            ;;
    esac
    echo "$phase complete: $phase_dir"
}

if [[ "${1:-}" == __worker ]]; then
    [[ "$#" -eq 4 ]] || die "internal worker argument error"
    run_worker "$2" "$4" "$3"
    exit 0
fi

if [[ "${1:-}" == mark-applied ]]; then
    [[ "$#" -eq 3 ]] || die "usage: $0 mark-applied PHASE OUT_DIR"
    mark_applied "$2" "$3"
    exit 0
fi

if [[ "$#" -ne 2 ]]; then
    echo "usage: $0 {dense-count|small|dense-seeds|peel|heavy} OUT_DIR" >&2
    exit 2
fi
run_phase "$1" "$2"
