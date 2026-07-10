#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

bash "$ROOT/experiments/peeling/build.sh" >/dev/null
bash "$ROOT/experiments/tiling/build.sh" >/dev/null

XOR="$ROOT/experiments/peeling/xor_bench"
TILING="$ROOT/experiments/tiling/rowop_tiling"

"$XOR" --sizes 1280 --target-gib 0.001 --repeats 3 \
    --working-mib 1 --seed 7 > "$work/default.csv" 2> "$work/default.err"
"$XOR" --cold --sizes 1280 --target-gib 0.001 --repeats 3 \
    --pool-gib 0.01 --seed 7 > "$work/cold.csv" 2> "$work/cold.err"
"$XOR" --muladd --sizes 1280 --target-gib 0.001 --repeats 3 \
    --pool-gib 0.01 --seed 7 > "$work/muladd.csv" 2> "$work/muladd.err"
"$XOR" --fanin 2 --sizes 1280 --target-gib 0.001 --repeats 3 \
    --pool-gib 0.01 --seed 7 > "$work/fanin.csv" 2> "$work/fanin.err"
if "$XOR" --cold --muladd --sizes 1280 --target-gib 0.001 --repeats 1 \
    --pool-gib 0.01 > "$work/mixed.csv" 2> "$work/mixed.err"; then
    echo "multiple CSV output modes unexpectedly succeeded" >&2
    exit 1
fi
if [[ -s "$work/mixed.csv" ]]; then
    echo "rejected output modes wrote a partial CSV" >&2
    exit 1
fi

cat > "$work/xor_schedule.csv" <<'EOF'
stage,op_type,dst_kind,dst_block,dst_offset,src0_kind,src0_block,src0_offset,src1_kind,src1_block,src1_offset,scalar,bytes
fixture,xor,recovery,0,0,input,0,0,none,-1,0,0,1280
EOF
"$TILING" --schedule "$work/xor_schedule.csv" --block-bytes 1280 \
    --tiles 256 --repeats 3 > "$work/xor_trace.csv" 2> "$work/xor_trace.err"
"$TILING" --schedule \
    "$ROOT/experiments/tiling/results/wh_oplog_N128_bb1280_full_20260703.csv" \
    --block-bytes 1280 --tiles 256 --repeats 3 \
    > "$work/real_trace.csv" 2> "$work/real_trace.err"

cat > "$work/overflow_schedule.csv" <<'EOF'
stage,op_type,dst_kind,dst_block,dst_offset,src0_kind,src0_block,src0_offset,src1_kind,src1_block,src1_offset,scalar,bytes
fixture,zero,recovery,1,0,none,-1,0,none,-1,0,0,1
EOF
overflow_block_bytes=$(python3 - <<'PY'
import struct
print(1 << (8 * struct.calcsize("P") - 1))
PY
)
if "$TILING" --verify-only --schedule "$work/overflow_schedule.csv" \
    --block-bytes "$overflow_block_bytes" --tiles 1 \
    > "$work/overflow.csv" 2> "$work/overflow.err"; then
    echo "overflowing trace storage unexpectedly succeeded" >&2
    exit 1
fi
if [[ -s "$work/overflow.csv" ]] ||
    ! grep -q 'schedule storage size overflows size_t' "$work/overflow.err"; then
    echo "overflowing trace storage did not fail cleanly before CSV output" >&2
    exit 1
fi

if "$TILING" --verify-only --schedule "$work/overflow_schedule.csv" \
    --block-bytes $((50 * 1024 * 1024)) --tiles 1 \
    > "$work/oversized.csv" 2> "$work/oversized.err"; then
    echo "oversized trace storage unexpectedly succeeded" >&2
    exit 1
fi
if [[ -s "$work/oversized.csv" ]] ||
    ! grep -q 'schedule aggregate storage exceeds --max-memory-mib policy' \
        "$work/oversized.err"; then
    echo "oversized trace storage did not fail cleanly before CSV output" >&2
    exit 1
fi

if "$TILING" --verify-only --schedule "$work/overflow_schedule.csv" \
    --block-bytes 1 --tiles 1 --max-memory-mib 1 \
    > "$work/op_storage.csv" 2> "$work/op_storage.err"; then
    echo "trace parser exceeded its aggregate memory budget" >&2
    exit 1
fi
if [[ -s "$work/op_storage.csv" ]] ||
    ! grep -q 'schedule operation storage exceeds --max-memory-mib policy' \
        "$work/op_storage.err"; then
    echo "trace operation storage policy did not fail cleanly" >&2
    exit 1
fi

cat > "$work/missing_header.csv" <<'EOF'
fixture,zero,recovery,0,0,none,-1,0,none,-1,0,0,1
fixture,zero,recovery,0,0,none,-1,0,none,-1,0,0,1
EOF
if "$TILING" --verify-only --schedule "$work/missing_header.csv" \
    --block-bytes 1 --tiles 1 \
    > "$work/missing_header.out" 2> "$work/missing_header.err"; then
    echo "headerless trace unexpectedly succeeded" >&2
    exit 1
fi
if [[ -s "$work/missing_header.out" ]] ||
    ! grep -q 'unexpected trace CSV header' "$work/missing_header.err"; then
    echo "headerless trace did not fail cleanly before CSV output" >&2
    exit 1
fi

cat > "$work/numeric_overflow.csv" <<'EOF'
stage,op_type,dst_kind,dst_block,dst_offset,src0_kind,src0_block,src0_offset,src1_kind,src1_block,src1_offset,scalar,bytes
fixture,zero,recovery,999999999999999999999999999999999999,0,none,-1,0,none,-1,0,0,1
EOF
if "$TILING" --verify-only --schedule "$work/numeric_overflow.csv" \
    --block-bytes 1 --tiles 1 \
    > "$work/numeric_overflow.out" 2> "$work/numeric_overflow.err"; then
    echo "overflowing trace integer unexpectedly succeeded" >&2
    exit 1
fi
if [[ -s "$work/numeric_overflow.out" ]] ||
    ! grep -q 'invalid trace row' "$work/numeric_overflow.err"; then
    echo "overflowing trace integer did not fail cleanly before CSV output" >&2
    exit 1
fi

python3 - "$work/long_line.csv" <<'PY'
import sys

with open(sys.argv[1], "w", encoding="ascii") as handle:
    handle.write("x" * 4097 + "\n")
PY
if "$TILING" --verify-only --schedule "$work/long_line.csv" \
    --block-bytes 1 --tiles 1 \
    > "$work/long_line.out" 2> "$work/long_line.err"; then
    echo "overlong trace line unexpectedly succeeded" >&2
    exit 1
fi
if [[ -s "$work/long_line.out" ]] ||
    ! grep -q 'schedule line exceeds 4096-byte limit' \
        "$work/long_line.err"; then
    echo "overlong trace line did not fail cleanly before CSV output" >&2
    exit 1
fi

cat > "$work/too_many_ops.csv" <<'EOF'
stage,op_type,dst_kind,dst_block,dst_offset,src0_kind,src0_block,src0_offset,src1_kind,src1_block,src1_offset,scalar,bytes
fixture,zero,recovery,0,0,none,-1,0,none,-1,0,0,1
fixture,zero,recovery,0,0,none,-1,0,none,-1,0,0,1
fixture,zero,recovery,0,0,none,-1,0,none,-1,0,0,1
EOF
if "$TILING" --verify-only --schedule "$work/too_many_ops.csv" \
    --block-bytes 1 --tiles 1 --max-operations 2 \
    > "$work/too_many_ops.out" 2> "$work/too_many_ops.err"; then
    echo "over-limit trace operation schedule unexpectedly succeeded" >&2
    exit 1
fi
if [[ -s "$work/too_many_ops.out" ]] ||
    ! grep -q 'schedule exceeds the 2-operation replay limit' \
        "$work/too_many_ops.err"; then
    echo "over-limit trace did not fail cleanly before CSV output" >&2
    exit 1
fi

cat > "$work/overlapping_dst_src.csv" <<'EOF'
stage,op_type,dst_kind,dst_block,dst_offset,src0_kind,src0_block,src0_offset,src1_kind,src1_block,src1_offset,scalar,bytes
fixture,memcpy,recovery,0,1,recovery,0,0,none,-1,0,0,16
EOF
if "$TILING" --verify-only --schedule "$work/overlapping_dst_src.csv" \
    --block-bytes 32 --tiles 8 \
    > "$work/overlapping_dst_src.out" 2> "$work/overlapping_dst_src.err"; then
    echo "overlapping destination/source trace unexpectedly succeeded" >&2
    exit 1
fi
if [[ -s "$work/overlapping_dst_src.out" ]] ||
    ! grep -q 'overlapping operands are unsupported for op_type=memcpy' \
        "$work/overlapping_dst_src.err"; then
    echo "overlapping destination/source trace did not fail cleanly" >&2
    exit 1
fi

cat > "$work/overlapping_sources.csv" <<'EOF'
stage,op_type,dst_kind,dst_block,dst_offset,src0_kind,src0_block,src0_offset,src1_kind,src1_block,src1_offset,scalar,bytes
fixture,addset,recovery,0,0,input,0,0,input,0,1,0,16
EOF
if "$TILING" --verify-only --schedule "$work/overlapping_sources.csv" \
    --block-bytes 32 --tiles 8 \
    > "$work/overlapping_sources.out" 2> "$work/overlapping_sources.err"; then
    echo "overlapping two-source trace unexpectedly succeeded" >&2
    exit 1
fi
if [[ -s "$work/overlapping_sources.out" ]] ||
    ! grep -q 'overlapping operands are unsupported for op_type=addset' \
        "$work/overlapping_sources.err"; then
    echo "overlapping two-source trace did not fail cleanly" >&2
    exit 1
fi

cat > "$work/adjacent_operands.csv" <<'EOF'
stage,op_type,dst_kind,dst_block,dst_offset,src0_kind,src0_block,src0_offset,src1_kind,src1_block,src1_offset,scalar,bytes
fixture,memcpy,recovery,0,16,recovery,0,0,none,-1,0,0,16
EOF
if ! "$TILING" --verify-only --schedule "$work/adjacent_operands.csv" \
    --block-bytes 32 --tiles 8 \
    > "$work/adjacent_operands.out" 2> "$work/adjacent_operands.err" ||
    ! grep -qx 'schedule verify: ok' "$work/adjacent_operands.out"; then
    echo "adjacent non-overlapping trace did not remain valid" >&2
    exit 1
fi

if (ulimit -v 98304) 2>/dev/null; then
    set +e
    (
        ulimit -v 98304
        "$TILING" --schedule "$work/overflow_schedule.csv" \
            --block-bytes 2 --tiles 1 --repeats 16000000 \
            --max-memory-mib 160
    ) > "$work/allocation_failure.out" 2> "$work/allocation_failure.err"
    allocation_rc=$?
    set -e
    if [[ $allocation_rc -ne 1 || -s "$work/allocation_failure.out" ]] ||
        ! grep -q 'schedule replay allocation failed' \
            "$work/allocation_failure.err"; then
        echo "trace allocation failure emitted partial output or was unstable" >&2
        exit 1
    fi
fi

python3 "$ROOT/experiments/validate_byte_metrics.py" \
    "$ROOT/experiments/peeling/results/cold_xor_calibration.csv" \
    "$ROOT/experiments/peeling/results/muladd_calibration.csv" \
    "$ROOT/experiments/peeling/results/fanin_gather.csv" \
    "$ROOT/experiments/tiling/results/rowop_tiling_1280_100k_1m.csv" \
    "$ROOT/experiments/tiling/results/rowop_tiling_real_oplog_20260703.csv" \
    "$work/default.csv" "$work/cold.csv" "$work/muladd.csv" \
    "$work/fanin.csv" "$work/xor_trace.csv" "$work/real_trace.csv" \
    > "$work/validation.out"

WORK="$work" REPO_ROOT="$ROOT" python3 - <<'PY'
import csv
import os

root = os.environ["WORK"]

def rows(name):
    with open(os.path.join(root, name), newline="") as handle:
        return list(csv.DictReader(handle))

for name in ("default.csv", "cold.csv"):
    row = rows(name)[0]
    assert (int(row["logical_bytes_per_op"]),
            int(row["estimated_read_bytes_per_op"]),
            int(row["estimated_write_bytes_per_op"])) == (1280, 2560, 1280)

muladd = rows("muladd.csv")[0]
assert (int(muladd["logical_bytes_per_op"]),
        int(muladd["estimated_read_bytes_per_op"]),
        int(muladd["estimated_write_bytes_per_op"])) == (1280, 2560, 1280)

fanin = rows("fanin.csv")[0]
assert (int(fanin["chained_estimated_read_bytes_per_op"]),
        int(fanin["chained_estimated_write_bytes_per_op"])) == (5120, 2560)
assert (int(fanin["gather_estimated_read_bytes_per_op"]),
        int(fanin["gather_estimated_write_bytes_per_op"])) == (3840, 1280)

xor_trace = rows("xor_trace.csv")
assert len(xor_trace) == 2
assert len({row["checksum"] for row in xor_trace}) == 1
for row in xor_trace:
    assert abs(float(row["estimated_read_gib"]) -
               2 * float(row["logical_work_gib"])) < 2e-6
    assert abs(float(row["estimated_write_gib"]) -
               float(row["logical_work_gib"])) < 2e-6

real_trace = rows("real_trace.csv")
assert len(real_trace) == 2
assert {row["checksum"] for row in real_trace} == {"0x8c08d158d46a7cbe"}
times = [float(row["median_ms"]) for row in real_trace]
assert all(value > 0.0 for value in times)
speed_ratio = times[0] / times[1]
with open(os.path.join(
        os.environ["REPO_ROOT"],
        "experiments/tiling/results/rowop_tiling_real_oplog_20260703.csv"),
        newline="") as handle:
    historical = list(csv.DictReader(handle))
old_rows = [row for row in historical
            if row["case"].endswith("oplog_N128_bb1280_full") and
            (row["mode"] == "untiled" or row["tile_bytes"] == "256")]
assert len(old_rows) == 2
old_rows.sort(key=lambda row: row["mode"] != "untiled")
historical_ratio = (float(old_rows[0]["median_ms"]) /
                    float(old_rows[1]["median_ms"]))
assert historical_ratio / 3.0 < speed_ratio < historical_ratio * 3.0
print("byte-metric E2E passed; real-trace median ratio=%.3f (v1 %.3f)" %
      (speed_ratio, historical_ratio))
PY
