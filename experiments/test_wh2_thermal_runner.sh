#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
# shellcheck source-path=SCRIPTDIR
source "$ROOT/experiments/wh2_thermal_runner.sh"

work=$(mktemp -d /tmp/wh2-thermal-runner-test.XXXXXX)
trap 'rm -rf "$work"' EXIT

fixture=$work/thermal.csv
valid_start='2026-01-01T00:00:00Z,100.25,100,4000,60,40,41,42,43,44,45,46,47,0,1,1,1,0,0'
missing='2026-01-01T00:00:01Z,101.25,100,4000,61,,,,,,,,,8,1,1,1,0,0'
valid_end='2026-01-01T00:00:02Z,102.25,99.5,3950,62,50,51,52,53,54,55,56,57,0,2,2,2,0,0'

# CRLF intentionally exercises the historical delimiter-loss path: assigning
# to $NF must not let awk rebuild the record with its default space OFS.
printf '%s\r\n%s\r\n%s\r\n' \
    "$WH2_THERMAL_CSV_HEADER" "$valid_start" "$missing" >"$fixture"

latest=$(latest_valid_thermal "$fixture" 101.25)
[ "$latest" = "$valid_start" ]
[ "$(awk -F, '{print NF}' <<<"$latest")" -eq 19 ]
case "$latest" in
    *,*,*) ;;
    *) printf '%s\n' 'latest thermal row lost CSV delimiters' >&2; exit 1 ;;
esac

start=$(thermal_monotonic_from_row "$latest")
wh2_thermal_numeric_scalar "$start"
[ "$start" = 100.25 ]
negative_mono=${valid_start/,100.25,/,\-0.25,}
if thermal_monotonic_from_row "$negative_mono" >/dev/null; then
    printf '%s\n' 'negative monotonic timestamp unexpectedly passed' >&2
    exit 1
fi

printf '%s\r\n' "$valid_end" >>"$fixture"
latest=$(latest_valid_thermal "$fixture" 102.25)
[ "$latest" = "$valid_end" ]
end=$(thermal_monotonic_from_row "$latest")
wh2_thermal_numeric_scalar "$end"
[ "$end" = 102.25 ]
awk -v start="$start" -v end="$end" 'BEGIN { exit !(start <= end) }'

window=$work/final/thermal_window.csv
summary=$work/final/thermal_summary.txt
mkdir -p "$work/final"
finalize_wh2_runner_telemetry 0 "$fixture" "$start" "$end" \
    "$window" "$summary"
[ -s "$window" ] && [ -s "$summary" ]
[ "$(awk 'END {print NR}' "$window")" -eq 4 ]
[ "$(awk -F, 'NR>1 && NF!=19 {++bad} END {print bad+0}' "$window")" -eq 0 ]
expected_window=$work/expected-window.csv
printf '%s\n%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" \
    "$missing" "$valid_end" >"$expected_window"
cmp -s "$expected_window" "$window"
[ "$(awk -F= '$1=="rows" {print $2}' "$summary")" = 3 ]
[ "$(awk -F= '$1=="valid_rows" {print $2}' "$summary")" = 2 ]
[ "$(awk -F= '$1=="missing_read_rows" {print $2}' "$summary")" = 1 ]
[ "$(awk -F= '$1=="invalid_rows" {print $2}' "$summary")" = 1 ]
[ "$(awk -F= '$1=="cpu_tctl_max_c" {print $2}' "$summary")" = 62 ]
[ "$(awk -F= '$1=="dimm_max_c" {print $2}' "$summary")" = 57 ]
if [ -n "$(find "$work/final" -maxdepth 1 -type f -name '*.part.*' \
    -print -quit)" ]; then
    printf '%s\n' 'thermal finalization leaked a staging file' >&2
    exit 1
fi

# Replacing a complete artifact pair must take the successful update path,
# publish the same validated contents, and remove both rollback copies.
replacement_window=$work/replacement-expected-window.csv
replacement_summary=$work/replacement-expected-summary.txt
cp -f "$window" "$replacement_window"
cp -f "$summary" "$replacement_summary"
printf '%s\n' stale-window >"$window"
printf '%s\n' stale-summary >"$summary"
finalize_wh2_runner_telemetry 0 "$fixture" "$start" "$end" \
    "$window" "$summary"
cmp -s "$replacement_window" "$window"
cmp -s "$replacement_summary" "$summary"
if [ -n "$(find "$work/final" -maxdepth 1 -type f \
    \( -name '*.part.*' -o -name '*.rollback.*' \) -print -quit)" ]; then
    printf '%s\n' 'successful replacement leaked a staging file' >&2
    exit 1
fi

# A sustained missing-DIMM interval is not fresh enough to stand in for a
# valid sample, while the transient row above was counted and excluded.
printf '%s\r\n' \
    '2026-01-01T00:00:08Z,108.25,100,4000,61,,,,,,,,,8,1,1,1,0,0' \
    >>"$fixture"
if latest_valid_thermal "$fixture" 108.25 >/dev/null; then
    printf '%s\n' 'sustained missing DIMM rows unexpectedly passed' >&2
    exit 1
fi
repeated_header=$work/repeated-header.csv
printf '%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" \
    "$WH2_THERMAL_CSV_HEADER" >"$repeated_header"
if latest_valid_thermal "$repeated_header" 100.25 >/dev/null; then
    printf '%s\n' 'repeated CSV header unexpectedly passed freshness' >&2
    exit 1
fi

# The freshness reader intentionally considers exactly the newest eight data
# rows.  A valid sample in the eighth slot must remain usable when followed by
# seven recent tolerated gaps; this distinguishes the lower boundary from a
# silently shortened tail window.
tail_lower=$work/tail-lower-boundary.csv
printf '%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" >"$tail_lower"
for mono in 100.75 101.25 101.75 102.25 102.75 103.25 103.75; do
    tail_missing_row=${missing/,101.25,/,$mono,}
    printf '%s\n' "$tail_missing_row" >>"$tail_lower"
done
[ "$(latest_valid_thermal "$tail_lower" 103.75)" = "$valid_start" ]

# Conversely, an older malformed row is irrelevant once it falls outside the
# eight-row window, while eight newer gaps push the last valid row out.
tail_window=$work/tail-window.csv
tail_old_invalid='2026-01-01T00:00:00Z,99,100,4000,oops,40,41,42,43,44,45,46,47,0,1,1,1,0,0'
printf '%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$tail_old_invalid" \
    >"$tail_window"
for mono in {101..108}; do
    tail_valid=${valid_start/,100.25,/,$mono,}
    printf '%s\n' "$tail_valid" >>"$tail_window"
done
latest=$(latest_valid_thermal "$tail_window" 108)
[ "$latest" = "$tail_valid" ]

tail_missing=$work/tail-missing.csv
printf '%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" \
    >"$tail_missing"
for mono in {101..108}; do
    tail_missing_row=${missing/,101.25,/,$mono,}
    printf '%s\n' "$tail_missing_row" >>"$tail_missing"
done
if latest_valid_thermal "$tail_missing" 108 >/dev/null; then
    printf '%s\n' 'valid sample outside tail window unexpectedly passed' >&2
    exit 1
fi

# Exercise the default CLOCK_BOOTTIME-derived freshness bound and the
# well-formed nonzero EDAC rejection path.
uptime_now=$(awk '{print $1}' /proc/uptime)
default_now_row=${valid_start/,100.25,/,$uptime_now,}
default_now_csv=$work/default-now.csv
printf '%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$default_now_row" \
    >"$default_now_csv"
[ "$(latest_valid_thermal "$default_now_csv")" = "$default_now_row" ]
# Put CE=1 on an otherwise tolerated missing-DIMM row after a prior valid row.
# It is rejected specifically by the sticky bad_edac latch; without that latch
# the freshness fallback would incorrectly return the prior valid row.
nonzero_edac_row=${missing/%,0,0/,1,0}
nonzero_edac_csv=$work/nonzero-edac-latest.csv
printf '%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" \
    "$nonzero_edac_row" >"$nonzero_edac_csv"
if latest_valid_thermal "$nonzero_edac_csv" 101.25 >/dev/null; then
    printf '%s\n' 'nonzero EDAC freshness sample unexpectedly passed' >&2
    exit 1
fi

# Header validation and data consumption must use one open file description.
# Rotate a live symlink from one internally valid snapshot to another exactly
# when the downstream reader starts; the result must still describe the file
# whose header was validated.
fd_original=$work/fd-original.csv
fd_replacement=$work/fd-replacement.csv
fd_live=$work/fd-live.csv
printf '%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" >"$fd_original"
printf '%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_end" >"$fd_replacement"
ln -s "$fd_original" "$fd_live"
fd_wrappers=$work/fd-wrappers
mkdir -p "$fd_wrappers"
fd_marker=$work/fd-wrapper-fired
real_tail=$(command -v tail)
printf '%s\n' '#!/bin/sh' \
    ": > \"\$FD_MARKER\"" \
    "ln -sfn \"\$FD_REPLACEMENT\" \"\$FD_LIVE\"" \
    "exec \"\$REAL_READER\" \"\$@\"" >"$fd_wrappers/tail"
chmod 755 "$fd_wrappers/tail"
fd_latest=$(FD_LIVE="$fd_live" FD_REPLACEMENT="$fd_replacement" \
    FD_MARKER="$fd_marker" REAL_READER="$real_tail" PATH="$fd_wrappers:$PATH" \
    latest_valid_thermal "$fd_live" 100.25)
[ "$fd_latest" = "$valid_start" ]
[ -e "$fd_marker" ]
rm -f "$fd_marker"
ln -sfn "$fd_original" "$fd_live"
real_awk=$(command -v awk)
cp -f "$fd_wrappers/tail" "$fd_wrappers/awk"
fd_summary=$(FD_LIVE="$fd_live" FD_REPLACEMENT="$fd_replacement" \
    FD_MARKER="$fd_marker" REAL_READER="$real_awk" PATH="$fd_wrappers:$PATH" \
    summarize_thermal_csv "$fd_live")
[ "$(awk -F= '$1=="cpu_tctl_max_c" {print $2}' \
    <<<"$fd_summary")" = 60 ]
[ -e "$fd_marker" ]
rm -f "$fd_marker"
rm -f "$fd_wrappers/awk"
printf '%s\n' '#!/bin/sh' \
    "\"\$REAL_READER\" \"\$@\"" 'exit 7' >"$fd_wrappers/tail"
chmod 755 "$fd_wrappers/tail"
for pipefail_mode in off on; do
    if (
        if [ "$pipefail_mode" = on ]; then
            set -o pipefail
        else
            set +o pipefail
        fi
        REAL_READER="$real_tail" PATH="$fd_wrappers:$PATH" \
            latest_valid_thermal "$fd_original" 100.25 >/dev/null
    ); then
        printf '%s\n' "failed tail passed with pipefail $pipefail_mode" >&2
        exit 1
    fi
done

# Codec-analysis failure and post-analysis telemetry failure must remain
# distinguishable even though both happen after all worker jobs have exited.
analysis_window=$work/analysis-window.csv
analysis_summary=$work/analysis-summary.txt
if finalize_wh2_runner_telemetry 1 "$fixture" "$start" "$end" \
    "$analysis_window" "$analysis_summary"
then
    printf '%s\n' 'codec-analysis failure unexpectedly succeeded' >&2
    exit 1
else
    rc=$?
fi
[ "$rc" -eq "$WH2_CODEC_ANALYSIS_FAILURE_RC" ]
[ ! -e "$analysis_window" ] && [ ! -e "$analysis_summary" ]

header_only=$work/header-only.csv
printf '%s\n' "$WH2_THERMAL_CSV_HEADER" >"$header_only"
telemetry_window=$work/telemetry-window.csv
telemetry_summary=$work/telemetry-summary.txt
if finalize_wh2_runner_telemetry 0 "$header_only" 1 2 \
    "$telemetry_window" "$telemetry_summary"
then
    printf '%s\n' 'empty telemetry window unexpectedly succeeded' >&2
    exit 1
else
    rc=$?
fi
[ "$rc" -eq "$WH2_POST_ANALYSIS_TELEMETRY_FAILURE_RC" ]
[ ! -e "$telemetry_window" ] && [ ! -e "$telemetry_summary" ]
if [ -n "$(find "$work" -maxdepth 1 -type f -name '*.part.*' \
    -print -quit)" ]; then
    printf '%s\n' 'failed thermal finalization leaked a staging file' >&2
    exit 1
fi

# Staging failures must not expose empty/partial replacements or disturb a
# previously committed artifact pair.
printf '%s\n' old-window >"$telemetry_window"
printf '%s\n' old-summary >"$telemetry_summary"
if finalize_wh2_runner_telemetry 0 "$header_only" 1 2 \
    "$telemetry_window" "$telemetry_summary"
then
    printf '%s\n' 'failed thermal finalization replaced existing artifacts' >&2
    exit 1
else
    rc=$?
fi
[ "$rc" -eq "$WH2_POST_ANALYSIS_TELEMETRY_FAILURE_RC" ]
[ "$(cat "$telemetry_window")" = old-window ]
[ "$(cat "$telemetry_summary")" = old-summary ]

# Lexical aliases and directory destinations must fail before publication.
alias_dir=$work/alias
mkdir -p "$alias_dir"
if build_thermal_artifacts "$fixture" "$start" "$end" \
    "$alias_dir/window" "$alias_dir/./window"
then
    printf '%s\n' 'lexically aliased output paths unexpectedly succeeded' >&2
    exit 1
fi
mkdir -p "$alias_dir/csv-dir" "$alias_dir/summary-dir"
if build_thermal_artifacts "$fixture" "$start" "$end" \
    "$alias_dir/csv-dir" "$alias_dir/summary"
then
    printf '%s\n' 'CSV directory destination unexpectedly succeeded' >&2
    exit 1
fi
if build_thermal_artifacts "$fixture" "$start" "$end" \
    "$alias_dir/window" "$alias_dir/summary-dir"
then
    printf '%s\n' 'summary directory destination unexpectedly succeeded' >&2
    exit 1
fi
fixture_copy=$work/fixture-copy.csv
cp -f "$fixture" "$fixture_copy"
if build_thermal_artifacts "$fixture" "$start" "$end" \
    "$fixture" "$alias_dir/source-summary"
then
    printf '%s\n' 'source/output alias unexpectedly succeeded' >&2
    exit 1
fi
cmp -s "$fixture_copy" "$fixture"
mkfifo "$alias_dir/source-fifo" "$alias_dir/output-fifo"
if build_thermal_artifacts "$alias_dir/source-fifo" "$start" "$end" \
    "$alias_dir/fifo-window" "$alias_dir/fifo-summary"
then
    printf '%s\n' 'FIFO telemetry source unexpectedly succeeded' >&2
    exit 1
fi
if build_thermal_artifacts "$fixture" "$start" "$end" \
    "$alias_dir/output-fifo" "$alias_dir/fifo-summary"
then
    printf '%s\n' 'FIFO output destination unexpectedly succeeded' >&2
    exit 1
fi
ln -s "$alias_dir/window-target" "$alias_dir/output-link"
if build_thermal_artifacts "$fixture" "$start" "$end" \
    "$alias_dir/output-link" "$alias_dir/link-summary"
then
    printf '%s\n' 'symlink output destination unexpectedly succeeded' >&2
    exit 1
fi

# A first-rename failure must leave the old CSV in place and restore its
# invalidated summary before returning failure.
rollback_window=$work/rollback-window.csv
rollback_summary=$work/rollback-summary.txt
printf '%s\n' old-window >"$rollback_window"
printf '%s\n' old-summary >"$rollback_summary"
mv_calls=0
mv()
{
    mv_calls=$((mv_calls + 1))
    if [ "$mv_calls" -eq 1 ]; then
        return 1
    fi
    command mv "$@"
}
if build_thermal_artifacts "$fixture" "$start" "$end" \
    "$rollback_window" "$rollback_summary"
then
    printf '%s\n' 'injected CSV-rename failure unexpectedly succeeded' >&2
    exit 1
fi
unset -f mv
[ "$(cat "$rollback_window")" = old-window ]
[ "$(cat "$rollback_summary")" = old-summary ]
if [ -n "$(find "$work" -maxdepth 1 -type f \
    \( -name '*.part.*' -o -name '*.rollback.*' \) -print -quit)" ]; then
    printf '%s\n' 'first-rename failure leaked a staging file' >&2
    exit 1
fi

# A second-rename failure must synchronously restore the previous pair.
mv_calls=0
mv()
{
    mv_calls=$((mv_calls + 1))
    if [ "$mv_calls" -eq 2 ]; then
        return 1
    fi
    command mv "$@"
}
if build_thermal_artifacts "$fixture" "$start" "$end" \
    "$rollback_window" "$rollback_summary"
then
    printf '%s\n' 'injected summary-rename failure unexpectedly succeeded' >&2
    exit 1
fi
unset -f mv
[ "$(cat "$rollback_window")" = old-window ]
[ "$(cat "$rollback_summary")" = old-summary ]
if [ -n "$(find "$work" -maxdepth 1 -type f \
    \( -name '*.part.*' -o -name '*.rollback.*' \) -print -quit)" ]; then
    printf '%s\n' 'rollback finalization leaked a staging file' >&2
    exit 1
fi

# If CSV rollback itself fails, the old summary stays invalidated and both
# recovery copies remain available rather than exposing a mismatched pair.
printf '%s\n' old-window >"$rollback_window"
printf '%s\n' old-summary >"$rollback_summary"
mv_calls=0
mv()
{
    mv_calls=$((mv_calls + 1))
    if [ "$mv_calls" -eq 2 ] || [ "$mv_calls" -eq 3 ]; then
        return 1
    fi
    command mv "$@"
}
if build_thermal_artifacts "$fixture" "$start" "$end" \
    "$rollback_window" "$rollback_summary"
then
    printf '%s\n' 'injected rollback failure unexpectedly succeeded' >&2
    exit 1
fi
unset -f mv
[ ! -e "$rollback_summary" ]
[ "$(find "$work" -maxdepth 1 -type f -name '*.rollback.*' | \
    wc -l)" -eq 2 ]
rm -f -- "$work"/*.rollback.* "$rollback_window"

# Sparse telemetry cannot claim continuous under-load coverage.
sparse=$work/sparse.csv
printf '%s\n' "$WH2_THERMAL_CSV_HEADER" >"$sparse"
printf '%s\n' "${valid_start/,100.25,/,100,}" >>"$sparse"
printf '%s\n' "${valid_end/,102.25,/,200,}" >>"$sparse"
if summarize_thermal_csv "$sparse" >/dev/null; then
    printf '%s\n' '100-second telemetry gap unexpectedly succeeded' >&2
    exit 1
fi

# Captured bounds must still be present as valid samples when the source is
# sliced; a rotated or truncated sampler file cannot silently shorten a run.
bounds_window=$work/bounds-window.csv
bounds_summary=$work/bounds-summary.txt
if build_thermal_artifacts "$fixture" 99 102.25 \
    "$bounds_window" "$bounds_summary"
then
    printf '%s\n' 'missing start-bound sample unexpectedly succeeded' >&2
    exit 1
fi

# EDAC counters are nonnegative integers, and loads are nonnegative scalars.
negative_edac=$work/negative-edac.csv
printf '%s\n' "$WH2_THERMAL_CSV_HEADER" >"$negative_edac"
printf '%s\n' "${valid_start/%,0,0/,-1,0}" >>"$negative_edac"
printf '%s\n' "$valid_end" >>"$negative_edac"
if summarize_thermal_csv "$negative_edac" >/dev/null; then
    printf '%s\n' 'negative EDAC counter unexpectedly succeeded' >&2
    exit 1
fi
negative_load=${missing/,1,1,1,0,0/,-1,1,1,0,0}
printf '%s\n%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" \
    "$negative_load" "$valid_end" \
    >"$work/negative-load.csv"
if summarize_thermal_csv "$work/negative-load.csv" >/dev/null; then
    printf '%s\n' 'negative load value unexpectedly succeeded' >&2
    exit 1
fi

# Only empty DIMM fields paired with a positive sampler read-error count are
# transient.  Nonnumeric values, zero-error blanks, and malformed CPU fields
# must not be hidden between otherwise-valid samples or by freshness fallback.
malformed_dimm='2026-01-01T00:00:01Z,101.25,100,4000,61,oops,,,,,,,,8,1,1,1,0,0'
zero_error_blank='2026-01-01T00:00:01Z,101.25,100,4000,61,,,,,,,,,0,1,1,1,0,0'
malformed_cpu='2026-01-01T00:00:01Z,101.25,100,4000,oops,40,41,42,43,44,45,46,47,0,1,1,1,0,0'
for row in "$malformed_dimm" "$zero_error_blank" "$malformed_cpu"; do
    malformed_csv=$work/malformed-$RANDOM.csv
    printf '%s\n%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" \
        "$valid_start" "$row" "$valid_end" >"$malformed_csv"
    if summarize_thermal_csv "$malformed_csv" >/dev/null; then
        printf '%s\n' 'malformed telemetry row unexpectedly summarized' >&2
        exit 1
    fi
    printf '%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" \
        "$valid_start" "$row" >"$malformed_csv"
    if latest_valid_thermal "$malformed_csv" 101.25 >/dev/null; then
        printf '%s\n' 'malformed latest row unexpectedly fell back' >&2
        exit 1
    fi
done
malformed_mono=${missing/,101.25,/,oops,}
printf '%s\n%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" \
    "$malformed_mono" "$valid_end" >"$work/malformed-mono-source.csv"
if build_thermal_artifacts "$work/malformed-mono-source.csv" \
    "$start" "$end" "$work/malformed-mono-window.csv" \
    "$work/malformed-mono-summary.txt"
then
    printf '%s\n' 'malformed in-window timestamp unexpectedly disappeared' >&2
    exit 1
fi
for bad_mono in 99 1000; do
    bad_missing=${missing/,101.25,/,${bad_mono},}
    printf '%s\n%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" \
        "$bad_missing" "$valid_end" >"$work/bad-missing-mono.csv"
    if build_thermal_artifacts "$work/bad-missing-mono.csv" \
        "$start" "$end" "$work/bad-missing-window.csv" \
        "$work/bad-missing-summary.txt"
    then
        printf '%s\n' 'out-of-order missing-DIMM timestamp passed' >&2
        exit 1
    fi
done
far_missing=${missing/,101.25,/,1000,}
printf '%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" \
    "$far_missing" >"$work/far-missing.csv"
if summarize_thermal_csv "$work/far-missing.csv" >/dev/null; then
    printf '%s\n' 'large missing-DIMM sampling gap unexpectedly passed' >&2
    exit 1
fi
# CLOCK_BOOTTIME samples must advance strictly.  Equal timestamps are neither
# ordered evidence nor a second independent thermal observation.
duplicate_mono=${valid_end/,102.25,/,101.25,}
printf '%s\n%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" \
    "$duplicate_mono" "$duplicate_mono" >"$work/duplicate-mono.csv"
if summarize_thermal_csv "$work/duplicate-mono.csv" >/dev/null; then
    printf '%s\n' 'duplicate monotonic timestamps unexpectedly passed' >&2
    exit 1
fi
fractional_errors=${missing/,8,1,1,1,0,0/,1.5,1,1,1,0,0}
printf '%s\n%s\n%s\n%s\n' "$WH2_THERMAL_CSV_HEADER" "$valid_start" \
    "$fractional_errors" "$valid_end" >"$work/fractional-errors.csv"
if summarize_thermal_csv "$work/fractional-errors.csv" >/dev/null; then
    printf '%s\n' 'fractional DIMM read-error count unexpectedly passed' >&2
    exit 1
fi

# Missing arguments remain classified even when the caller uses nounset.
if ( set -u; finalize_wh2_runner_telemetry ); then
    printf '%s\n' 'missing codec RC unexpectedly succeeded' >&2
    exit 1
else
    rc=$?
fi
[ "$rc" -eq "$WH2_CODEC_ANALYSIS_FAILURE_RC" ]
if ( set -u; finalize_wh2_runner_telemetry 0 ); then
    printf '%s\n' 'missing artifact arguments unexpectedly succeeded' >&2
    exit 1
else
    rc=$?
fi
[ "$rc" -eq "$WH2_POST_ANALYSIS_TELEMETRY_FAILURE_RC" ]
if finalize_wh2_runner_telemetry \
    999999999999999999999999999999999999999999999999 \
    "$fixture" "$start" "$end" "$window" "$summary"
then
    printf '%s\n' 'oversized nonzero codec RC unexpectedly succeeded' >&2
    exit 1
else
    rc=$?
fi
[ "$rc" -eq "$WH2_CODEC_ANALYSIS_FAILURE_RC" ]

printf '%s\n' \
    'wh2 thermal runner: CSV/scalars/window/missing-DIMM/failure-classification PASS'
