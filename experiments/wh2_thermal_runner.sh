#!/usr/bin/env bash

# Thermal-window helpers shared by long-running WH2 benchmark runners.
# Source this file; the functions deliberately do not change the caller's
# shell options.

WH2_THERMAL_CSV_HEADER='utc,monotonic_s,cpu_busy_pct,cpu_avg_mhz,cpu_tctl_c,dimm_i2c1_50_c,dimm_i2c1_51_c,dimm_i2c1_52_c,dimm_i2c1_53_c,dimm_i2c2_50_c,dimm_i2c2_51_c,dimm_i2c2_52_c,dimm_i2c2_53_c,dimm_read_errors,load1,load5,load15,edac_ce,edac_ue'
WH2_CODEC_ANALYSIS_FAILURE_RC=89
WH2_POST_ANALYSIS_TELEMETRY_FAILURE_RC=90

wh2_thermal_numeric_scalar()
{
    [[ ${1:-} =~ ^[0-9]+([.][0-9]+)?$ ]]
}

latest_valid_thermal()
{
    [ "$#" -ge 1 ] && [ "$#" -le 2 ] || return 1
    local path=$1
    local now=${2:-}
    local header thermal_fd
    local -a pipeline_status

    if [ -z "$now" ]; then
        now=$(awk '{print $1}' /proc/uptime) || return 1
    fi
    wh2_thermal_numeric_scalar "$now" || return 1
    [ -s "$path" ] || return 1
    exec {thermal_fd}<"$path" || return 1
    IFS= read -r header <&"$thermal_fd" || {
        exec {thermal_fd}<&-
        return 1
    }
    header=${header%$'\r'}
    if [ "$header" != "$WH2_THERMAL_CSV_HEADER" ]; then
        exec {thermal_fd}<&-
        return 1
    fi

    tail -n 8 <&"$thermal_fd" | \
        LC_ALL=C awk -F, -v now="$now" '
        BEGIN { OFS=FS }
        function number(v) { return v ~ /^-?[0-9]+([.][0-9]+)?$/ }
        function counter(v) { return v ~ /^[0-9]+$/ }
        {
            # Assigning to a field makes awk rebuild $0.  OFS must remain a
            # comma or a CRLF input row silently becomes space-separated.
            sub(/\r$/, "", $NF)
            if (NF != 19 || !number($2)) {
                hard_invalid=1
                next
            }
            if (!have_newest || $2 > newest) newest=$2
            have_newest=1
            other_valid=$2>=0 && number($3) && $3>=0 && $3<=100 &&
                number($4) && $4>0 && number($5) && $5>=0 && $5<=120 &&
                counter($14) &&
                number($15) && $15>=0 && number($16) && $16>=0 &&
                number($17) && $17>=0 && counter($18) && counter($19)
            dimm_missing=0
            dimm_bad=0
            for (i=6; i<=13; ++i) {
                if ($i == "") dimm_missing=1
                else if (!number($i) || $i<0 || $i>100) dimm_bad=1
            }
            tolerated_missing=other_valid && $14>0 &&
                dimm_missing && !dimm_bad
            valid=other_valid && $14==0 && !dimm_missing && !dimm_bad &&
                $18==0 && $19==0
            if (counter($18) && counter($19) &&
                ($18 != 0 || $19 != 0)) bad_edac=1
            if (!valid && !tolerated_missing) hard_invalid=1
            if (valid && (!have_valid || $2 >= valid_mono)) {
                valid_mono=$2
                valid_line=$0
                have_valid=1
            }
        }
        END {
            if (bad_edac || hard_invalid || !have_newest || !have_valid ||
                newest-valid_mono > 5.0 || now-valid_mono > 5.0 ||
                valid_mono-now > 1.0)
                exit 1
            print valid_line
        }
    '
    pipeline_status=("${PIPESTATUS[@]}")
    exec {thermal_fd}<&-
    if [ "${pipeline_status[0]}" -ne 0 ] ||
       [ "${pipeline_status[1]}" -ne 0 ]
    then
        return 1
    fi
}

thermal_monotonic_from_row()
{
    [ "$#" -eq 1 ] || return 1
    local row=$1

    printf '%s\n' "$row" | LC_ALL=C awk -F, '
        function number(v) { return v ~ /^[0-9]+([.][0-9]+)?$/ }
        {
            sub(/\r$/, "", $NF)
            if (NR != 1 || NF != 19 || !number($2) || $2 < 0) exit 1
            value=$2
        }
        END {
            if (NR != 1 || value == "") exit 1
            print value
        }
    '
}

summarize_thermal_csv()
{
    [ "$#" -eq 1 ] || [ "$#" -eq 3 ] || return 1
    local path=$1
    local expected_start=${2:-}
    local expected_end=${3:-}
    local header thermal_fd rc

    if [ "$#" -eq 3 ]; then
        wh2_thermal_numeric_scalar "$expected_start" || return 1
        wh2_thermal_numeric_scalar "$expected_end" || return 1
    fi

    [ -s "$path" ] || return 1
    exec {thermal_fd}<"$path" || return 1
    IFS= read -r header <&"$thermal_fd" || {
        exec {thermal_fd}<&-
        return 1
    }
    header=${header%$'\r'}
    if [ "$header" != "$WH2_THERMAL_CSV_HEADER" ]; then
        exec {thermal_fd}<&-
        return 1
    fi

    LC_ALL=C awk -F, -v expected_start="$expected_start" \
        -v expected_end="$expected_end" '
        BEGIN { OFS=FS }
        function number(v) { return v ~ /^-?[0-9]+([.][0-9]+)?$/ }
        function counter(v) { return v ~ /^[0-9]+$/ }
        {
            sub(/\r$/, "", $NF)
            ++rows
            if (NF==19 && number($2)) {
                if (have_row_mono) {
                    if ($2 <= row_mono) row_order_invalid=1
                    if ($2-row_mono > 5.0) row_sampling_gap=1
                }
                if (expected_start != "" &&
                    ($2 < expected_start || $2 > expected_end))
                    row_bounds_invalid=1
                row_mono=$2
                have_row_mono=1
            }
            other_valid=NF==19 && number($2) && $2>=0 &&
                number($3) && $3>=0 && $3<=100 &&
                number($4) && $4>0 &&
                number($5) && $5>=0 && $5<=120 &&
                counter($14) &&
                number($15) && $15>=0 && number($16) && $16>=0 &&
                number($17) && $17>=0 &&
                counter($18) && counter($19)
            dimm_missing=0
            dimm_bad=0
            for (i=6; i<=13; ++i) {
                if (NF==19 && $i == "") dimm_missing=1
                else if (NF!=19 || !number($i) || $i<0 || $i>100)
                    dimm_bad=1
            }
            missing=NF==19 &&
                ((counter($14) && $14!=0) || dimm_missing)
            if (missing) ++missing_rows
            if (NF!=19 || !counter($18) || !counter($19)) {
                valid=0
                edac_invalid=1
            } else {
                if (!have_ce || $18>ce_max) ce_max=$18
                if (!have_ue || $19>ue_max) ue_max=$19
                have_ce=have_ue=1
            }
            tolerated_missing=other_valid && $14>0 &&
                dimm_missing && !dimm_bad
            valid=other_valid && $14==0 && !dimm_missing && !dimm_bad
            if (!valid) {
                ++invalid_rows
                # Transient DIMM read gaps are the one intentionally
                # tolerated invalid-row class.  Any other malformed or
                # out-of-range telemetry must make the artifact unusable.
                if (!tolerated_missing) hard_invalid=1
                next
            }
            ++valid_rows
            if (valid_rows == 1) first_valid_mono=$2
            if (have_valid_mono &&
                ($2 <= valid_mono || $2-valid_mono > 5.0))
                sampling_gap=1
            valid_mono=$2
            have_valid_mono=1
            if (valid_rows==1 || $5>cpu_max) cpu_max=$5
            for (i=6; i<=13; ++i) {
                if (!have_dimm || $i>dimm_max) dimm_max=$i
                have_dimm=1
            }
        }
        END {
            print "schema=wirehair.wh2.thermal.v1"
            print "rows=" rows+0
            print "valid_rows=" valid_rows+0
            print "missing_read_rows=" missing_rows+0
            print "invalid_rows=" invalid_rows+0
            print "cpu_tctl_max_c=" cpu_max+0
            print "dimm_max_c=" dimm_max+0
            print "edac_ce_max=" ce_max+0
            print "edac_ue_max=" ue_max+0
            if (rows==0 || valid_rows==0 || sampling_gap ||
                row_sampling_gap || hard_invalid ||
                row_order_invalid || row_bounds_invalid || edac_invalid ||
                (expected_start != "" &&
                    (first_valid_mono != expected_start ||
                     valid_mono != expected_end)) ||
                !have_ce || !have_ue || ce_max!=0 || ue_max!=0)
                exit 1
        }
    ' <&"$thermal_fd"
    rc=$?
    exec {thermal_fd}<&-
    return "$rc"
}

build_thermal_artifacts()
{
    [ "$#" -eq 5 ] || return 1
    local source_csv=$1
    local start=$2
    local end=$3
    local output_csv=$4
    local output_summary=$5
    local source_real output_csv_real output_summary_real
    local tmp_csv tmp_summary rollback_csv='' rollback_summary=''
    local csv_restored

    wh2_thermal_numeric_scalar "$start" || return 1
    wh2_thermal_numeric_scalar "$end" || return 1
    [ -f "$source_csv" ] || return 1
    [ ! -d "$output_csv" ] && [ ! -d "$output_summary" ] || return 1
    if [ -L "$output_csv" ] ||
       { [ -e "$output_csv" ] && [ ! -f "$output_csv" ]; } ||
       [ -L "$output_summary" ] ||
       { [ -e "$output_summary" ] && [ ! -f "$output_summary" ]; }
    then
        return 1
    fi
    source_real=$(realpath -e -- "$source_csv") || return 1
    output_csv_real=$(realpath -m -- "$output_csv") || return 1
    output_summary_real=$(realpath -m -- "$output_summary") || return 1
    [ "$source_real" != "$output_csv_real" ] &&
        [ "$source_real" != "$output_summary_real" ] &&
        [ "$output_csv_real" != "$output_summary_real" ] || return 1
    if { [ -e "$output_csv" ] || [ -L "$output_csv" ]; } &&
       [ "$source_csv" -ef "$output_csv" ]
    then
        return 1
    fi
    if { [ -e "$output_summary" ] || [ -L "$output_summary" ]; } &&
       [ "$source_csv" -ef "$output_summary" ]
    then
        return 1
    fi
    if { [ -e "$output_csv" ] || [ -L "$output_csv" ]; } &&
       { [ -e "$output_summary" ] || [ -L "$output_summary" ]; } &&
       [ "$output_csv" -ef "$output_summary" ]
    then
        return 1
    fi
    LC_ALL=C awk -v start="$start" -v end="$end" \
        'BEGIN { exit !(start <= end) }' || return 1

    tmp_csv=$(mktemp "${output_csv}.part.XXXXXX") || return 1
    tmp_summary=$(mktemp "${output_summary}.part.XXXXXX") || {
        rm -f -- "$tmp_csv"
        return 1
    }

    if ! LC_ALL=C awk -F, -v start="$start" -v end="$end" '
        BEGIN { OFS=FS }
        {
            sub(/\r$/, "", $NF)
            if (NR==1) { print; next }
            numeric=$2 ~ /^[0-9]+([.][0-9]+)?$/
            if (!inside) {
                if (!numeric || $2 != start) next
                inside=1
                print
                if ($2 == end) { found_end=1; exit }
                next
            }
            # Once the exact start sample is found, retain every intervening
            # record.  Malformed timestamps must reach the validator rather
            # than disappearing from an apparently continuous window.
            print
            if (numeric && $2 == end) { found_end=1; exit }
        }
        END { if (!inside || !found_end) exit 1 }
    ' "$source_real" >"$tmp_csv" ||
       ! summarize_thermal_csv \
            "$tmp_csv" "$start" "$end" >"$tmp_summary"
    then
        rm -f -- "$tmp_csv" "$tmp_summary"
        return 1
    fi

    # Validate both staged files before publishing either.  Keep a synchronous
    # rollback copy of an existing window because the summary is the final
    # commit record and its rename can still fail after the CSV rename.
    if [ -e "$output_csv" ] || [ -L "$output_csv" ]
    then
        rollback_csv=$(mktemp "${output_csv}.rollback.XXXXXX") || {
            rm -f -- "$tmp_csv" "$tmp_summary"
            return 1
        }
        if ! cp -pf -- "$output_csv" "$rollback_csv"; then
            rm -f -- "$tmp_csv" "$tmp_summary" "$rollback_csv"
            return 1
        fi
    fi
    if [ -e "$output_summary" ] || [ -L "$output_summary" ]
    then
        rollback_summary=$(mktemp \
            "${output_summary}.rollback.XXXXXX") || {
            rm -f -- "$tmp_csv" "$tmp_summary" "$rollback_csv"
            return 1
        }
        if ! cp -pf -- "$output_summary" "$rollback_summary"; then
            rm -f -- "$tmp_csv" "$tmp_summary" \
                "$rollback_csv" "$rollback_summary"
            return 1
        fi
        # The summary is the commit record.  Invalidate it before replacing
        # the CSV so readers never accept a new window with an old summary.
        if ! rm -f -- "$output_summary"; then
            rm -f -- "$tmp_csv" "$tmp_summary" \
                "$rollback_csv" "$rollback_summary"
            return 1
        fi
    fi
    if ! mv -fT -- "$tmp_csv" "$output_csv"; then
        rm -f -- "$tmp_csv" "$tmp_summary" "$rollback_csv"
        if [ -n "$rollback_summary" ] &&
           ! mv -fT -- "$rollback_summary" "$output_summary"
        then
            # The prior CSV is still in place, but leave the recovery copy
            # when restoring its commit record is impossible.
            return 1
        fi
        return 1
    fi
    if ! mv -fT -- "$tmp_summary" "$output_summary"; then
        rm -f -- "$tmp_summary"
        csv_restored=1
        if [ -n "$rollback_csv" ]; then
            if ! mv -fT -- "$rollback_csv" "$output_csv"; then
                csv_restored=0
            fi
        elif ! rm -f -- "$output_csv"; then
            csv_restored=0
        fi
        # Restore the old summary only after the old CSV is back.  A failed
        # rollback therefore leaves no stale commit record and retains its
        # .rollback files for explicit recovery.
        if [ "$csv_restored" -eq 1 ] &&
           [ -n "$rollback_summary" ] &&
           ! mv -fT -- "$rollback_summary" "$output_summary"
        then
            return 1
        fi
        if [ "$csv_restored" -eq 0 ]; then
            return 1
        fi
        return 1
    fi
    rm -f -- "$rollback_csv" "$rollback_summary"
}

finalize_wh2_runner_telemetry()
{
    if [ "$#" -lt 1 ]; then
        return "$WH2_CODEC_ANALYSIS_FAILURE_RC"
    fi
    local codec_analysis_rc=$1
    shift

    if ! [[ "$codec_analysis_rc" =~ ^0+$ ]]; then
        return "$WH2_CODEC_ANALYSIS_FAILURE_RC"
    fi
    if [ "$#" -ne 5 ]; then
        return "$WH2_POST_ANALYSIS_TELEMETRY_FAILURE_RC"
    fi
    build_thermal_artifacts "$@" ||
        return "$WH2_POST_ANALYSIS_TELEMETRY_FAILURE_RC"
}
