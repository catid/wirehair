#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build"}
OUT_DIR=${OUT_DIR:-"$ROOT_DIR/tables/results/slice_$(date -u +%Y%m%dT%H%M%SZ)"}
REPEAT=${REPEAT:-1}

if [[ "$BUILD_DIR" != /* ]]; then
    BUILD_DIR="$ROOT_DIR/$BUILD_DIR"
fi
if [[ "$OUT_DIR" != /* ]]; then
    OUT_DIR="$ROOT_DIR/$OUT_DIR"
fi

SEED=${SEED:-123456789}

DCOUNT_TRIALS=${DCOUNT_TRIALS:-8}
DCOUNT_NLO=${DCOUNT_NLO:-2}
DCOUNT_NHI=${DCOUNT_NHI:-16}

PEEL_TRIALS=${PEEL_TRIALS:-1}
PEEL_SUBLO=${PEEL_SUBLO:-0}
PEEL_SUBHI=${PEEL_SUBHI:-0}
PEEL_NLO=${PEEL_NLO:-2048}
PEEL_NHI=${PEEL_NHI:-2048}
PEEL_MAX_TRIES=${PEEL_MAX_TRIES:-1}
PEEL_SKIP_TUNING=${PEEL_SKIP_TUNING:-0}

MOST_DENSE_TRIALS=${MOST_DENSE_TRIALS:-1}
MOST_DENSE_INDEX_LO=${MOST_DENSE_INDEX_LO:-13}
MOST_DENSE_INDEX_HI=${MOST_DENSE_INDEX_HI:-13}

SMALL_DENSE_TRIALS=${SMALL_DENSE_TRIALS:-1}
SMALL_DENSE_NLO=${SMALL_DENSE_NLO:-17}
SMALL_DENSE_NHI=${SMALL_DENSE_NHI:-17}

GEN_TABLES_NO_BENCHMARKS=${GEN_TABLES_NO_BENCHMARKS:-1}
HEAVY_TRIALS=${HEAVY_TRIALS:-0}

if [[ -z "${OMP_NUM_THREADS:-}" ]]; then
    OMP_NUM_THREADS=$(nproc 2>/dev/null || echo 1)
    export OMP_NUM_THREADS
fi

mkdir -p "$OUT_DIR"
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

    if [[ "$REPEAT" == 1 ]]; then
        "$@" > "$OUT_DIR/$name.repeat.out" 2> "$OUT_DIR/$name.repeat.err"
        diff -u "$OUT_DIR/$name.out" "$OUT_DIR/$name.repeat.out" \
            > "$OUT_DIR/$name.repeat.out.diff"
        diff -u "$OUT_DIR/$name.err" "$OUT_DIR/$name.repeat.err" \
            > "$OUT_DIR/$name.repeat.err.diff"
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
    echo "omp_num_threads=$OMP_NUM_THREADS"
    echo "seed=$SEED"
    echo
} > "$MANIFEST"

build_generators

run_step dense_count \
    "$BUILD_DIR/gen_dcounts" \
    --seed "$SEED" \
    --trials "$DCOUNT_TRIALS" \
    --nlo "$DCOUNT_NLO" \
    --nhi "$DCOUNT_NHI"

peel_cmd=(
    "$BUILD_DIR/gen_peel_seeds"
    --trials "$PEEL_TRIALS"
    --sublo "$PEEL_SUBLO"
    --subhi "$PEEL_SUBHI"
    --nlo "$PEEL_NLO"
    --nhi "$PEEL_NHI"
    --max-tries "$PEEL_MAX_TRIES"
)
if [[ "$PEEL_SKIP_TUNING" == 1 ]]; then
    peel_cmd+=(--skip-tuning)
fi
run_step peel_seeds "${peel_cmd[@]}"

run_step most_dense_seeds \
    "$BUILD_DIR/gen_most_dseeds" \
    --seed "$SEED" \
    --trials "$MOST_DENSE_TRIALS" \
    --dense-index-lo "$MOST_DENSE_INDEX_LO" \
    --dense-index-hi "$MOST_DENSE_INDEX_HI"

run_step small_dense_seeds \
    "$BUILD_DIR/gen_small_dseeds" \
    --seed "$SEED" \
    --trials "$SMALL_DENSE_TRIALS" \
    --nlo "$SMALL_DENSE_NLO" \
    --nhi "$SMALL_DENSE_NHI"

gen_tables_cmd=("$BUILD_DIR/gen_tables" --heavy-trials "$HEAVY_TRIALS")
if [[ "$GEN_TABLES_NO_BENCHMARKS" == 1 ]]; then
    gen_tables_cmd+=(--no-benchmarks)
fi
run_step gen_tables "${gen_tables_cmd[@]}"

echo "slice complete: $OUT_DIR"
