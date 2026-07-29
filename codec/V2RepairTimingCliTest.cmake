if(NOT DEFINED BENCH)
    message(FATAL_ERROR "BENCH is required")
endif()

set(clean_env
    "${CMAKE_COMMAND}" -E env
    --unset=WIREHAIR_V2_PEEL_DEGREES
    --unset=WIREHAIR_V2_STAIRCASE_DEGREES
    --unset=WIREHAIR_V2_STAIRCASE_ROW_DEGREES
    --unset=WIREHAIR_V2_STAIRCASE_DEGREE_SCALE
    --unset=WIREHAIR_V2_BAND_TRACKING_X)

set(base_args
    repairtiming
    --N 16
    --bb 2
    --repair-arm pure8_s0m1_d3_repair_v1
    --repair-policy repair-v1
    --dispatch-profile dispatch-v1
    --construction-seed-derivation derived
    --construction-seed 123
    --loss 0.1
    --loss-seed 456
    --schedule adversarial
    --warmup-replicates 0
    --replicates 3
    --inner-reps 1
    --max-overhead 64
    --cache-state warm
    --systematic-cache off
    --evict-bytes 4096
    --context-sha256
    0000000000000000000000000000000000000000000000000000000000000000
    --required-margin 0)

execute_process(
    COMMAND ${clean_env} "${BENCH}" ${base_args}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err
    TIMEOUT 45)
if(NOT result EQUAL 0 OR NOT err STREQUAL "")
    message(FATAL_ERROR
        "repairtiming base command failed\n"
        "rc=${result}\nstdout=${out}\nstderr=${err}")
endif()

string(CONCAT done_regex
    "\n# repairtiming_done,complete=1,cells=3,summary_rows=3,"
    "attempt_rows=24,timing_rows=336,rows=363,"
    "finished_monotonic_ns=[0-9]+,stream_sha256=[0-9a-f]+\n$")
if(NOT out MATCHES
       "^# repairtiming,protocol=wirehair-v2-bench:repairtiming:repair-v1:v3,"
   OR NOT out MATCHES
       ",timing_stats_policy=first-inner-structural-counters-summed-phase-ns-v1,selector_timing_stats_policy=all-executed-selector-subcalls-summed-per-inner-v1,structural_digest_encoding="
   OR NOT out MATCHES
       "\n# repairtiming_summary_columns,row_kind,cell_index,"
   OR NOT out MATCHES
       "\n# repairtiming_attempt_columns,row_kind,cell_index,"
   OR NOT out MATCHES
       "\n# repairtiming_timing_columns,row_kind,cell_index,"
   OR NOT out MATCHES "${done_regex}")
    message(FATAL_ERROR
        "repairtiming envelope/cardinality mismatch\n${out}")
endif()

string(REGEX MATCH
    "# repairtiming_done,[^\n]*\n$" done "${out}")
string(REGEX REPLACE
    "# repairtiming_done,[^\n]*\n$" "" body "${out}")
string(REGEX REPLACE
    "stream_sha256=[0-9a-f]+\n$" "stream_sha256="
    done_prefix "${done}")
string(SHA256 bound_sha "${body}${done_prefix}")
string(REGEX MATCH
    "stream_sha256=([0-9a-f]+)" hash_match "${done}")
set(done_sha "${CMAKE_MATCH_1}")
if(NOT bound_sha STREQUAL done_sha)
    message(FATAL_ERROR
        "repairtiming stream SHA-256 mismatch")
endif()

# The v3 domain is the full finite half-open interval [0, 1), rather than
# bandtiming's older [0, 0.99] convenience cap.
set(high_loss_args ${base_args})
list(FIND high_loss_args "0.1" high_loss_index)
if(high_loss_index LESS 0)
    message(FATAL_ERROR "repairtiming test lost its --loss value")
endif()
list(REMOVE_AT high_loss_args ${high_loss_index})
list(INSERT high_loss_args ${high_loss_index} "0.995")
execute_process(
    COMMAND ${clean_env} "${BENCH}" ${high_loss_args}
    RESULT_VARIABLE high_loss_result
    OUTPUT_VARIABLE high_loss_out
    ERROR_VARIABLE high_loss_err
    TIMEOUT 45)
if(NOT high_loss_result EQUAL 0 OR
   NOT high_loss_err STREQUAL "" OR
   NOT high_loss_out MATCHES
       "^# repairtiming,[^\n]*,loss=0.995[0-9]*,")
    message(FATAL_ERROR
        "repairtiming half-open loss domain failed\n"
        "rc=${high_loss_result}\nstdout=${high_loss_out}\n"
        "stderr=${high_loss_err}")
endif()

string(REGEX MATCH
    "# repairtiming_summary_columns,([^\n]+)"
    summary_header_match "${out}")
set(summary_header "${CMAKE_MATCH_1}")
string(REPLACE "," ";" summary_fields "${summary_header}")
list(LENGTH summary_fields summary_width)
string(REGEX MATCH
    "# repairtiming_attempt_columns,([^\n]+)"
    attempt_header_match "${out}")
set(attempt_header "${CMAKE_MATCH_1}")
string(REPLACE "," ";" attempt_fields "${attempt_header}")
list(LENGTH attempt_fields attempt_width)
string(REGEX MATCH
    "# repairtiming_timing_columns,([^\n]+)"
    timing_header_match "${out}")
set(timing_header "${CMAKE_MATCH_1}")
string(REPLACE "," ";" timing_fields "${timing_header}")
list(LENGTH timing_fields timing_width)
if(NOT summary_width EQUAL 80 OR
   NOT attempt_width EQUAL 86 OR
   NOT timing_width EQUAL 71)
    message(FATAL_ERROR
        "repairtiming typed widths changed: "
        "${summary_width}/${attempt_width}/${timing_width}")
endif()

string(REGEX MATCHALL "\nsummary,[^\n]+" summary_rows "${out}")
string(REGEX MATCHALL "\nattempt,[^\n]+" attempt_rows "${out}")
string(REGEX MATCHALL "\ntiming,[^\n]+" timing_rows "${out}")
list(LENGTH summary_rows summary_count)
list(LENGTH attempt_rows attempt_count)
list(LENGTH timing_rows timing_count)
if(NOT summary_count EQUAL 3 OR
   NOT attempt_count EQUAL 24 OR
   NOT timing_count EQUAL 336)
    message(FATAL_ERROR
        "repairtiming row-kind counts changed: "
        "${summary_count}/${attempt_count}/${timing_count}")
endif()

list(FIND summary_fields forced_equal forced_equal_index)
list(FIND summary_fields repair_direct_executed repair_direct_index)
list(FIND summary_fields dispatch_direct_executed dispatch_direct_index)
list(FIND summary_fields direct_fixed_prefix_symbols prefix_index)
list(FIND summary_fields selected_attempt selected_attempt_index)
list(FIND summary_fields descriptor_hex descriptor_hex_index)
list(FIND summary_fields descriptor_sha256 descriptor_sha_index)
list(FIND summary_fields repaired_recovery_ok repaired_ok_index)
list(FIND summary_fields repair_direct_result repair_direct_result_index)
list(FIND summary_fields
    repair_direct_witness_equal repair_direct_witness_index)
list(FIND summary_fields
    dispatch_intermediate_bytes dispatch_intermediate_index)
list(FIND summary_fields calls_executed calls_executed_index)
list(FIND summary_fields
    real_configuration_calls real_calls_index)
list(FIND summary_fields
    structural_probe_calls probe_calls_index)
list(GET summary_rows 0 first_summary)
string(SUBSTRING "${first_summary}" 1 -1 first_summary)
string(REPLACE "," ";" first_summary_fields "${first_summary}")
list(GET first_summary_fields ${forced_equal_index} forced_equal)
list(GET first_summary_fields ${repair_direct_index} repair_direct_executed)
list(GET first_summary_fields ${dispatch_direct_index}
    dispatch_direct_executed)
list(GET first_summary_fields ${prefix_index} direct_prefix)
if(NOT forced_equal STREQUAL "1" OR
   NOT repair_direct_executed STREQUAL "1" OR
   NOT dispatch_direct_executed STREQUAL "1" OR
   direct_prefix LESS 16 OR direct_prefix GREATER 80)
    message(FATAL_ERROR
        "repairtiming semantic bridge/direct prefix failed\n"
        "${first_summary}")
endif()

# Construction root zero is a frozen, payload-independent weak pure8 case at
# K=16: attempt zero is singular and repair-v1 commits attempt one.  Exercise
# the actual serialized repair path rather than testing healthy attempt zero
# alone.
set(repaired_args ${base_args})
list(FIND repaired_args
    "--construction-seed-derivation" repaired_derivation_option)
list(FIND repaired_args "--construction-seed" repaired_seed_option)
if(repaired_derivation_option LESS 0 OR repaired_seed_option LESS 0)
    message(FATAL_ERROR "repairtiming test lost construction options")
endif()
math(EXPR repaired_derivation_value
    "${repaired_derivation_option} + 1")
math(EXPR repaired_seed_value "${repaired_seed_option} + 1")
list(REMOVE_AT repaired_args ${repaired_derivation_value})
list(INSERT repaired_args ${repaired_derivation_value} "fixed")
list(REMOVE_AT repaired_args ${repaired_seed_value})
list(INSERT repaired_args ${repaired_seed_value} "0")
execute_process(
    COMMAND ${clean_env} "${BENCH}" ${repaired_args}
    RESULT_VARIABLE repaired_result
    OUTPUT_VARIABLE repaired_out
    ERROR_VARIABLE repaired_err
    TIMEOUT 45)
if(NOT repaired_result EQUAL 0 OR
   NOT repaired_err STREQUAL "" OR
   NOT repaired_out MATCHES
       "construction_seed_derivation=fixed_base_u32_v1,")
    message(FATAL_ERROR
        "repairtiming repaired-root command failed\n"
        "rc=${repaired_result}\nstdout=${repaired_out}\n"
        "stderr=${repaired_err}")
endif()
string(REGEX MATCHALL "\nsummary,[^\n]+" repaired_summaries
    "${repaired_out}")
list(LENGTH repaired_summaries repaired_summary_count)
if(NOT repaired_summary_count EQUAL 3)
    message(FATAL_ERROR
        "repairtiming repaired-root summary count changed")
endif()
set(expected_descriptor_hex
    "5732523309bfe05c77cfcc190000000001")
foreach(repaired_summary IN LISTS repaired_summaries)
    string(SUBSTRING "${repaired_summary}" 1 -1 repaired_summary)
    string(REPLACE "," ";" repaired_values "${repaired_summary}")
    list(GET repaired_values
        ${selected_attempt_index} repaired_attempt)
    list(GET repaired_values ${descriptor_hex_index} repaired_descriptor)
    list(GET repaired_values
        ${descriptor_sha_index} repaired_descriptor_sha)
    list(GET repaired_values ${forced_equal_index} repaired_forced_equal)
    list(GET repaired_values ${repaired_ok_index} repaired_recovery_ok)
    list(GET repaired_values
        ${repair_direct_index} repaired_direct_executed)
    list(GET repaired_values
        ${repair_direct_result_index} repaired_direct_result)
    list(GET repaired_values
        ${repair_direct_witness_index} repaired_direct_witness)
    list(GET repaired_values
        ${dispatch_intermediate_index} repaired_dispatch_intermediate)
    list(GET repaired_values
        ${calls_executed_index} repaired_calls)
    list(GET repaired_values ${real_calls_index} repaired_real_calls)
    list(GET repaired_values ${probe_calls_index} repaired_probe_calls)
    if(NOT repaired_attempt STREQUAL "1" OR
       NOT repaired_descriptor STREQUAL "${expected_descriptor_hex}" OR
       NOT repaired_descriptor_sha MATCHES "^[0-9a-f]+$" OR
       NOT repaired_calls STREQUAL "4" OR
       NOT repaired_real_calls STREQUAL "2" OR
       NOT repaired_probe_calls STREQUAL "2" OR
       NOT repaired_forced_equal STREQUAL "1" OR
       NOT repaired_recovery_ok STREQUAL "1" OR
       NOT repaired_direct_executed STREQUAL "1" OR
       NOT repaired_direct_result STREQUAL "0" OR
       NOT repaired_direct_witness STREQUAL "1" OR
       NOT repaired_dispatch_intermediate STREQUAL "74")
        message(FATAL_ERROR
            "repairtiming repaired-root bridge failed\n"
            "${repaired_summary}")
    endif()
endforeach()

list(FIND attempt_fields probe_executed probe_executed_index)
list(FIND attempt_fields probe_block_xors probe_block_xors_index)
list(FIND attempt_fields real_executed real_executed_index)
list(FIND attempt_fields real_block_xors real_block_xors_index)
set(repaired_selector_block_xors_0 0)
set(repaired_selector_block_xors_1 0)
set(repaired_selector_block_xors_2 0)
string(REGEX MATCHALL "\nattempt,[^\n]+" repaired_attempts
    "${repaired_out}")
foreach(repaired_attempt_row IN LISTS repaired_attempts)
    string(SUBSTRING "${repaired_attempt_row}" 1 -1
        repaired_attempt_row)
    string(REPLACE "," ";" repaired_attempt_values
        "${repaired_attempt_row}")
    list(GET repaired_attempt_values 1 repaired_attempt_cell)
    list(GET repaired_attempt_values
        ${probe_executed_index} repaired_probe_executed)
    list(GET repaired_attempt_values
        ${real_executed_index} repaired_real_executed)
    set(repaired_call_block_xors 0)
    if(repaired_probe_executed STREQUAL "1")
        list(GET repaired_attempt_values
            ${probe_block_xors_index} repaired_probe_block_xors)
        math(EXPR repaired_call_block_xors
            "${repaired_call_block_xors} + ${repaired_probe_block_xors}")
    endif()
    if(repaired_real_executed STREQUAL "1")
        list(GET repaired_attempt_values
            ${real_block_xors_index} repaired_real_block_xors)
        math(EXPR repaired_call_block_xors
            "${repaired_call_block_xors} + ${repaired_real_block_xors}")
    endif()
    set(repaired_sum_var
        "repaired_selector_block_xors_${repaired_attempt_cell}")
    math(EXPR repaired_sum_value
        "${${repaired_sum_var}} + ${repaired_call_block_xors}")
    set(${repaired_sum_var} ${repaired_sum_value})
endforeach()

list(FIND timing_fields timing_panel timing_panel_index)
list(FIND timing_fields timing_role timing_role_index)
list(FIND timing_fields timing_scope timing_scope_index)
list(FIND timing_fields timing_executed timing_executed_index)
list(FIND timing_fields timing_block_xors timing_block_xors_index)
string(REGEX MATCHALL "\ntiming,[^\n]+" repaired_timings
    "${repaired_out}")
set(repaired_selector_timing_count 0)
set(repaired_wh1_encoder_timing_count 0)
foreach(repaired_timing IN LISTS repaired_timings)
    string(SUBSTRING "${repaired_timing}" 1 -1 repaired_timing)
    string(REPLACE "," ";" repaired_timing_values
        "${repaired_timing}")
    list(GET repaired_timing_values
        ${timing_panel_index} repaired_panel)
    list(GET repaired_timing_values
        ${timing_role_index} repaired_role)
    if(repaired_role STREQUAL "wh1_encoder")
        list(GET repaired_timing_values
            ${timing_scope_index} repaired_scope)
        if(NOT repaired_scope STREQUAL "encoder_wh1")
            message(FATAL_ERROR
                "repairtiming WH1 encoder scope changed: "
                "${repaired_scope}")
        endif()
        math(EXPR repaired_wh1_encoder_timing_count
            "${repaired_wh1_encoder_timing_count} + 1")
    endif()
    if(repaired_panel STREQUAL "encoder_selector_forced" AND
       repaired_role STREQUAL "repair_selector_encoder")
        list(GET repaired_timing_values
            ${timing_executed_index} repaired_timing_executed)
        if(NOT repaired_timing_executed STREQUAL "1")
            message(FATAL_ERROR
                "repairtiming repaired selector timing was not executed")
        endif()
        list(GET repaired_timing_values 1 repaired_timing_cell)
        list(GET repaired_timing_values
            ${timing_block_xors_index} repaired_timing_block_xors)
        set(repaired_sum_var
            "repaired_selector_block_xors_${repaired_timing_cell}")
        if(NOT repaired_timing_block_xors EQUAL
            ${${repaired_sum_var}})
            message(FATAL_ERROR
                "repairtiming selector counters omit an executed call: "
                "cell=${repaired_timing_cell} "
                "timed=${repaired_timing_block_xors} "
                "calls=${${repaired_sum_var}}")
        endif()
        math(EXPR repaired_selector_timing_count
            "${repaired_selector_timing_count} + 1")
    endif()
endforeach()
if(NOT repaired_selector_timing_count EQUAL 12)
    message(FATAL_ERROR
        "repairtiming repaired selector timing count changed: "
        "${repaired_selector_timing_count}")
endif()
if(NOT repaired_wh1_encoder_timing_count EQUAL 36)
    message(FATAL_ERROR
        "repairtiming WH1 encoder timing count changed: "
        "${repaired_wh1_encoder_timing_count}")
endif()

function(expect_no_output_failure label)
    execute_process(
        COMMAND ${clean_env} "${BENCH}" ${ARGN}
        RESULT_VARIABLE failure_result
        OUTPUT_VARIABLE failure_out
        ERROR_VARIABLE failure_err
        TIMEOUT 45)
    if(failure_result EQUAL 0 OR NOT failure_out STREQUAL "" OR
       failure_err STREQUAL "")
        message(FATAL_ERROR
            "${label}: expected fail-closed rejection\n"
            "rc=${failure_result}\nstdout=${failure_out}\n"
            "stderr=${failure_err}")
    endif()
endfunction()

foreach(fault IN ITEMS bad-magic bad-id bad-attempt trailing)
    expect_no_output_failure(
        "transactional descriptor ${fault}"
        ${base_args} --test-descriptor-fault ${fault})
endforeach()
expect_no_output_failure(
    "timed hash scope guard"
    ${base_args} --test-timed-scope-drift hash)
expect_no_output_failure(
    "forced intermediate witness guard"
    ${base_args} --test-forced-witness-drift intermediate)
set(nonfinite_margin_args ${base_args})
list(FIND nonfinite_margin_args "--required-margin" margin_option_index)
if(margin_option_index LESS 0)
    message(FATAL_ERROR "repairtiming test lost --required-margin")
endif()
math(EXPR margin_value_index "${margin_option_index} + 1")
list(REMOVE_AT nonfinite_margin_args ${margin_value_index})
list(INSERT nonfinite_margin_args ${margin_value_index} "nan")
expect_no_output_failure(
    "non-finite required margin"
    ${nonfinite_margin_args})
expect_no_output_failure(
    "non-frozen max overhead"
    repairtiming
    --N 16 --bb 2
    --repair-arm pure8_s0m1_d3_repair_v1
    --repair-policy repair-v1
    --dispatch-profile dispatch-v1
    --construction-seed-derivation derived
    --construction-seed 123 --loss 0.1 --loss-seed 456
    --schedule adversarial --warmup-replicates 0 --replicates 3
    --inner-reps 1 --max-overhead 63 --cache-state warm
    --systematic-cache off --evict-bytes 4096
    --context-sha256
    0000000000000000000000000000000000000000000000000000000000000000
    --required-margin 0)

foreach(width IN ITEMS 6 32 64 256 1280 4096)
    execute_process(
        COMMAND ${clean_env} "${BENCH}"
            repairtiming
            --N 17 --bb ${width}
            --repair-arm pure9_s0m1_d3_repair_v1
            --repair-policy repair-v1
            --dispatch-profile dispatch-v1
            --construction-seed-derivation fixed
            --construction-seed 789
            --loss 0.2 --loss-seed 987
            --schedule burst
            --warmup-replicates 0 --replicates 3
            --inner-reps 1 --max-overhead 64
            --cache-state warm --systematic-cache on
            --evict-bytes 4096
            --context-sha256
            1111111111111111111111111111111111111111111111111111111111111111
            --required-margin 0
        RESULT_VARIABLE width_result
        OUTPUT_VARIABLE width_out
        ERROR_VARIABLE width_err
        TIMEOUT 45)
    if(NOT width_result EQUAL 0 OR NOT width_err STREQUAL "" OR
       NOT width_out MATCHES
           "bb=${width},message_bytes=[0-9]+,message_tail_bytes=1,"
       OR NOT width_out MATCHES
           "construction_seed_derivation=fixed_base_u32_v1,"
       OR NOT width_out MATCHES
           "\n# repairtiming_done,complete=1,cells=3,")
        message(FATAL_ERROR
            "repairtiming width/tail/fixed-root gate failed bb=${width}\n"
            "rc=${width_result}\nstdout=${width_out}\nstderr=${width_err}")
    endif()
endforeach()

execute_process(
    COMMAND ${clean_env}
        WIREHAIR_V2_PEEL_DEGREES=1,1
        "${BENCH}" ${base_args}
    RESULT_VARIABLE env_result
    OUTPUT_VARIABLE env_out
    ERROR_VARIABLE env_err
    TIMEOUT 10)
if(env_result EQUAL 0 OR NOT env_out STREQUAL "" OR
   NOT env_err MATCHES
       "repairtiming forbids ambient WIREHAIR_V2_PEEL_DEGREES")
    message(FATAL_ERROR
        "repairtiming ambient-hook rejection failed\n"
        "rc=${env_result}\nstdout=${env_out}\nstderr=${env_err}")
endif()
