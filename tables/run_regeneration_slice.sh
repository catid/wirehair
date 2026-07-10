#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build"}
OUT_DIR=${OUT_DIR:-"$ROOT_DIR/tables/results/slice_$(date -u +%Y%m%dT%H%M%SZ)"}
REPEAT=${REPEAT:-1}
BUILD_GENERATORS=${BUILD_GENERATORS:-1}

RUN_DENSE_COUNT=${RUN_DENSE_COUNT:-1}
RUN_PEEL_SEEDS=${RUN_PEEL_SEEDS:-1}
RUN_MOST_DENSE_SEEDS=${RUN_MOST_DENSE_SEEDS:-1}
RUN_SMALL_DENSE_SEEDS=${RUN_SMALL_DENSE_SEEDS:-1}
RUN_GEN_TABLES=${RUN_GEN_TABLES:-1}

if [[ "$BUILD_DIR" != /* ]]; then
    BUILD_DIR="$ROOT_DIR/$BUILD_DIR"
fi
if [[ "$OUT_DIR" != /* ]]; then
    OUT_DIR="$ROOT_DIR/$OUT_DIR"
fi

SEED=${SEED:-123456789}

DCOUNT_SEED=${DCOUNT_SEED:-$SEED}
DCOUNT_TRIALS=${DCOUNT_TRIALS:-8}
DCOUNT_NLO=${DCOUNT_NLO:-2}
DCOUNT_NHI=${DCOUNT_NHI:-16}
DCOUNT_MAX_FAILURES=${DCOUNT_MAX_FAILURES:-10}
DCOUNT_LOW_COUNT_RUN=${DCOUNT_LOW_COUNT_RUN:-4}

PEEL_TRIALS=${PEEL_TRIALS:-1}
PEEL_SEED=${PEEL_SEED:-$SEED}
PEEL_SUBLO=${PEEL_SUBLO:-0}
PEEL_SUBHI=${PEEL_SUBHI:-0}
PEEL_NLO=${PEEL_NLO:-2048}
PEEL_NHI=${PEEL_NHI:-2048}
PEEL_MAX_TRIES=${PEEL_MAX_TRIES:-1}
PEEL_LOSS10_PERCENT=${PEEL_LOSS10_PERCENT:-10}
PEEL_LOSS30_PERCENT=${PEEL_LOSS30_PERCENT:-30}
PEEL_SKIP_TUNING=${PEEL_SKIP_TUNING:-0}

MOST_DENSE_TRIALS=${MOST_DENSE_TRIALS:-1}
MOST_DENSE_SEED=${MOST_DENSE_SEED:-$SEED}
MOST_DENSE_INDEX_LO=${MOST_DENSE_INDEX_LO:-13}
MOST_DENSE_INDEX_HI=${MOST_DENSE_INDEX_HI:-13}
MOST_DENSE_N_SAMPLES=${MOST_DENSE_N_SAMPLES:-2}
MOST_DENSE_FULL_SEARCH=${MOST_DENSE_FULL_SEARCH:-0}

SMALL_DENSE_TRIALS=${SMALL_DENSE_TRIALS:-1}
SMALL_DENSE_SEED=${SMALL_DENSE_SEED:-$SEED}
SMALL_DENSE_NLO=${SMALL_DENSE_NLO:-17}
SMALL_DENSE_NHI=${SMALL_DENSE_NHI:-17}
SMALL_DENSE_NLIST=${SMALL_DENSE_NLIST:-}
SMALL_DENSE_FULL_SEARCH=${SMALL_DENSE_FULL_SEARCH:-0}

GEN_TABLES_NO_BENCHMARKS=${GEN_TABLES_NO_BENCHMARKS:-1}
HEAVY_TRIALS=${HEAVY_TRIALS:-0}

if [[ -z "${OMP_NUM_THREADS:-}" ]]; then
    OMP_NUM_THREADS=$(nproc 2>/dev/null || echo 1)
    export OMP_NUM_THREADS
fi

if [[ -e "$OUT_DIR" ]]; then
    if [[ ! -d "$OUT_DIR" || -n "$(find "$OUT_DIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
        echo "OUT_DIR already exists and is not empty: $OUT_DIR" >&2
        exit 2
    fi
else
    mkdir -p "$OUT_DIR"
fi
MANIFEST="$OUT_DIR/manifest.txt"

build_generators() {
    local args=(
        --build "$BUILD_DIR"
        --target
        gen_dcounts
        gen_peel_seeds
        gen_most_dseeds
        gen_small_dseeds
        gen_tables
    )
    if [[ -n "${BUILD_JOBS:-}" ]]; then
        args+=(-- -j "$BUILD_JOBS")
    fi
    cmake "${args[@]}"
}

quote_command() {
    printf '%q ' "$@"
    printf '\n'
}

checksum_outputs() {
    local name=$1
    sha256sum "$OUT_DIR/$name.out" "$OUT_DIR/$name.err" >> "$MANIFEST"
}

run_step() {
    local name=$1
    shift

    quote_command "$@" > "$OUT_DIR/$name.cmd"
    "$@" > "$OUT_DIR/$name.out" 2> "$OUT_DIR/$name.err"
    checksum_outputs "$name"
    if [[ -f "$OUT_DIR/$name.results.tsv" ]]; then
        sha256sum "$OUT_DIR/$name.results.tsv" >> "$MANIFEST"
    fi

    if [[ "$REPEAT" == 1 ]]; then
        local repeat_args=("$@")
        local i
        for ((i = 0; i + 1 < ${#repeat_args[@]}; ++i)); do
            if [[ "${repeat_args[i]}" == --results ]]; then
                repeat_args[i + 1]="$OUT_DIR/$name.repeat.results.tsv"
            fi
        done
        "${repeat_args[@]}" > "$OUT_DIR/$name.repeat.out" \
            2> "$OUT_DIR/$name.repeat.err"
        diff -u "$OUT_DIR/$name.out" "$OUT_DIR/$name.repeat.out" \
            > "$OUT_DIR/$name.repeat.out.diff"
        diff -u "$OUT_DIR/$name.err" "$OUT_DIR/$name.repeat.err" \
            > "$OUT_DIR/$name.repeat.err.diff"
        if [[ -f "$OUT_DIR/$name.results.tsv" ]]; then
            diff -u "$OUT_DIR/$name.results.tsv" \
                "$OUT_DIR/$name.repeat.results.tsv" \
                > "$OUT_DIR/$name.repeat.results.tsv.diff"
            sha256sum "$OUT_DIR/$name.repeat.results.tsv" >> "$MANIFEST"
        fi
        sha256sum "$OUT_DIR/$name.repeat.out" "$OUT_DIR/$name.repeat.err" \
            >> "$MANIFEST"
    fi
}

{
    echo "# wirehair table regeneration slice"
    echo "utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "root=$ROOT_DIR"
    echo "build_dir=$BUILD_DIR"
    echo "out_dir=$OUT_DIR"
    echo "git_commit=$(git -C "$ROOT_DIR" rev-parse HEAD)"
    echo "repeat=$REPEAT"
    echo "build_generators=$BUILD_GENERATORS"
    echo "run_dense_count=$RUN_DENSE_COUNT"
    echo "run_peel_seeds=$RUN_PEEL_SEEDS"
    echo "run_most_dense_seeds=$RUN_MOST_DENSE_SEEDS"
    echo "run_small_dense_seeds=$RUN_SMALL_DENSE_SEEDS"
    echo "run_gen_tables=$RUN_GEN_TABLES"
    echo "omp_num_threads=$OMP_NUM_THREADS"
    echo "seed=$SEED"
    echo "dcount_seed=$DCOUNT_SEED"
    echo "dcount_trials=$DCOUNT_TRIALS"
    echo "dcount_nlo=$DCOUNT_NLO"
    echo "dcount_nhi=$DCOUNT_NHI"
    echo "dcount_max_failures=$DCOUNT_MAX_FAILURES"
    echo "dcount_low_count_run=$DCOUNT_LOW_COUNT_RUN"
    echo "dcount_selection=threshold-run-v1"
    echo "peel_trials=$PEEL_TRIALS"
    echo "peel_seed=$PEEL_SEED"
    echo "peel_sublo=$PEEL_SUBLO"
    echo "peel_subhi=$PEEL_SUBHI"
    echo "peel_nlo=$PEEL_NLO"
    echo "peel_nhi=$PEEL_NHI"
    echo "peel_max_tries=$PEEL_MAX_TRIES"
    echo "peel_loss10_percent=$PEEL_LOSS10_PERCENT"
    echo "peel_loss30_percent=$PEEL_LOSS30_PERCENT"
    echo "peel_skip_tuning=$PEEL_SKIP_TUNING"
    echo "most_dense_trials=$MOST_DENSE_TRIALS"
    echo "most_dense_seed=$MOST_DENSE_SEED"
    echo "most_dense_index_lo=$MOST_DENSE_INDEX_LO"
    echo "most_dense_index_hi=$MOST_DENSE_INDEX_HI"
    echo "most_dense_n_samples=$MOST_DENSE_N_SAMPLES"
    echo "most_dense_full_search=$MOST_DENSE_FULL_SEARCH"
    echo "small_dense_trials=$SMALL_DENSE_TRIALS"
    echo "small_dense_seed=$SMALL_DENSE_SEED"
    echo "small_dense_nlo=$SMALL_DENSE_NLO"
    echo "small_dense_nhi=$SMALL_DENSE_NHI"
    echo "small_dense_nlist=$SMALL_DENSE_NLIST"
    echo "small_dense_full_search=$SMALL_DENSE_FULL_SEARCH"
    echo "gen_tables_no_benchmarks=$GEN_TABLES_NO_BENCHMARKS"
    echo "heavy_trials=$HEAVY_TRIALS"
    echo
} > "$MANIFEST"

if [[ "$BUILD_GENERATORS" == 1 ]]; then
    build_generators
fi

if [[ "$RUN_DENSE_COUNT" == 1 ]]; then
    run_step dense_count \
        "$BUILD_DIR/gen_dcounts" \
        --seed "$DCOUNT_SEED" \
        --trials "$DCOUNT_TRIALS" \
        --nlo "$DCOUNT_NLO" \
        --nhi "$DCOUNT_NHI" \
        --max-failures "$DCOUNT_MAX_FAILURES" \
        --low-count-run "$DCOUNT_LOW_COUNT_RUN" \
        --results "$OUT_DIR/dense_count.results.tsv"
fi

if [[ "$RUN_PEEL_SEEDS" == 1 ]]; then
    peel_cmd=(
        "$BUILD_DIR/gen_peel_seeds"
        --seed "$PEEL_SEED"
        --trials "$PEEL_TRIALS"
        --sublo "$PEEL_SUBLO"
        --subhi "$PEEL_SUBHI"
        --nlo "$PEEL_NLO"
        --nhi "$PEEL_NHI"
        --max-tries "$PEEL_MAX_TRIES"
        --loss10-percent "$PEEL_LOSS10_PERCENT"
        --loss30-percent "$PEEL_LOSS30_PERCENT"
        --results "$OUT_DIR/peel_seeds.results.tsv"
    )
    if [[ "$PEEL_SKIP_TUNING" == 1 ]]; then
        peel_cmd+=(--skip-tuning)
    fi
    run_step peel_seeds "${peel_cmd[@]}"
fi

if [[ "$RUN_MOST_DENSE_SEEDS" == 1 ]]; then
    most_dense_cmd=(
        "$BUILD_DIR/gen_most_dseeds" \
        --seed "$MOST_DENSE_SEED" \
        --trials "$MOST_DENSE_TRIALS" \
        --n-samples "$MOST_DENSE_N_SAMPLES" \
        --results "$OUT_DIR/most_dense_seeds.results.tsv"
    )
    if [[ "$MOST_DENSE_FULL_SEARCH" == 1 ]]; then
        most_dense_cmd+=(--full-search)
    else
        most_dense_cmd+=(
            --dense-index-lo "$MOST_DENSE_INDEX_LO"
            --dense-index-hi "$MOST_DENSE_INDEX_HI"
        )
    fi
    run_step most_dense_seeds "${most_dense_cmd[@]}"
fi

if [[ "$RUN_SMALL_DENSE_SEEDS" == 1 ]]; then
    small_dense_cmd=(
        "$BUILD_DIR/gen_small_dseeds"
        --seed "$SMALL_DENSE_SEED"
        --trials "$SMALL_DENSE_TRIALS"
        --nlo "$SMALL_DENSE_NLO"
        --nhi "$SMALL_DENSE_NHI"
        --results "$OUT_DIR/small_dense_seeds.results.tsv"
    )
    if [[ "$SMALL_DENSE_FULL_SEARCH" == 1 ]]; then
        small_dense_cmd+=(--full-search)
    elif [[ -n "$SMALL_DENSE_NLIST" ]]; then
        small_dense_cmd+=(--nlist "$SMALL_DENSE_NLIST")
    fi
    run_step small_dense_seeds "${small_dense_cmd[@]}"
fi

if [[ "$RUN_GEN_TABLES" == 1 ]]; then
    gen_tables_cmd=("$BUILD_DIR/gen_tables" --heavy-trials "$HEAVY_TRIALS")
    if [[ "$GEN_TABLES_NO_BENCHMARKS" == 1 ]]; then
        gen_tables_cmd+=(--no-benchmarks)
    fi
    run_step gen_tables "${gen_tables_cmd[@]}"
fi

echo "slice complete: $OUT_DIR"
