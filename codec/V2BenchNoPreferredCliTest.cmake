if(NOT DEFINED BENCH)
    message(FATAL_ERROR "BENCH is required")
endif()
if(NOT DEFINED PEELTIMING_SUPPORTED)
    if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
        set(PEELTIMING_SUPPORTED 1)
    else()
        set(PEELTIMING_SUPPORTED 0)
    endif()
endif()

# Keep command-gating checks independent of any experiment hook inherited from
# the developer shell.  This target intentionally retains the test hooks; the
# build_policy_e2e test owns the real strict BUILD_TESTS=OFF compile/run gate.
set(clean_hook_env_command
    "${CMAKE_COMMAND}" -E env
    --unset=WIREHAIR_V2_PEEL_DEGREES
    --unset=WIREHAIR_V2_STAIRCASE_DEGREES
    --unset=WIREHAIR_V2_STAIRCASE_ROW_DEGREES
    --unset=WIREHAIR_V2_STAIRCASE_DEGREE_SCALE
    --unset=WIREHAIR_V2_BAND_TRACKING_X)

execute_process(
    COMMAND ${clean_hook_env_command} "${BENCH}"
    RESULT_VARIABLE usage_result
    OUTPUT_VARIABLE usage_out
    ERROR_VARIABLE usage_err
    TIMEOUT 10)
if(NOT usage_result EQUAL 1 OR
   NOT usage_err MATCHES "^usage: wirehair_v2_bench " OR
   usage_err MATCHES "preferredattempt|preferredtiming" OR
   NOT usage_err MATCHES "peeltiming")
    message(FATAL_ERROR
        "preferredattempt leaked into preferred-suppressed usage\n"
        "rc=${usage_result}\nstdout=${usage_out}\nstderr=${usage_err}")
endif()

# peeltiming is a peel-PMF protocol, not a preferred-attempt experiment.  It
# must retain the same complete native receipt when only preferred routing is
# compiled out.
if(PEELTIMING_SUPPORTED)
execute_process(
    COMMAND ${clean_hook_env_command} "${BENCH}" peeltiming
    --N 16 --bb 2 --target-profile dispatch-v1 --seed-policy raw
    --construction-seed 7 --loss 0.1 --loss-seed 9 --schedule iid
    --candidate-pmf 0,1 --candidate-scale identity
    --warmup-replicates 0 --replicates 4 --inner-reps 1
    --max-overhead 8 --cache-state warm --evict-bytes 4096
    --context-sha256
    0000000000000000000000000000000000000000000000000000000000000000
    --required-margin 0.05
    RESULT_VARIABLE peeltiming_result
    OUTPUT_VARIABLE peeltiming_out
    ERROR_VARIABLE peeltiming_err
    TIMEOUT 10)
string(REGEX MATCHALL
    "\n[0-3],1,(candidate_control|identity_aa),[0-7],"
    peeltiming_rows "${peeltiming_out}")
list(LENGTH peeltiming_rows peeltiming_row_count)
string(REGEX MATCH
    "# peeltiming_done,[^\n]*\n$" peeltiming_done "${peeltiming_out}")
string(REGEX REPLACE
    "# peeltiming_done,[^\n]*\n$" "" peeltiming_body "${peeltiming_out}")
string(REGEX REPLACE
    "stream_sha256=[0-9a-f]+\n$" "stream_sha256="
    peeltiming_done_prefix "${peeltiming_done}")
string(SHA256 peeltiming_bound_sha256
    "${peeltiming_body}${peeltiming_done_prefix}")
string(REGEX MATCH
    "stream_sha256=([0-9a-f]+)" peeltiming_hash_match
    "${peeltiming_done}")
set(peeltiming_done_sha256 "${CMAKE_MATCH_1}")
if(NOT peeltiming_result EQUAL 0 OR peeltiming_err OR
   NOT peeltiming_row_count EQUAL 64 OR
   NOT peeltiming_out MATCHES
       "^# peeltiming,schema=wirehair.wh2.peeltiming.v2," OR
   NOT peeltiming_out MATCHES
       "slot_prewarm=validated-plus-conditioning-matching-role-solves-same-cpu-before-cache-v1,.*stream_hash_scope=body-plus-done-prefix-v1," OR
   NOT peeltiming_out MATCHES
       "\n# peel_semantic,timed=0,.*full_recovery_equal=1,pass=1\n" OR
   NOT peeltiming_out MATCHES
       "\n# peeltiming_done,complete=1,rows=64,.*\n$" OR
   NOT peeltiming_bound_sha256 STREQUAL peeltiming_done_sha256)
    message(FATAL_ERROR
        "peeltiming disappeared or became incomplete in the "
        "preferred-suppressed build\nrc=${peeltiming_result}\n"
        "stdout=${peeltiming_out}\nstderr=${peeltiming_err}")
endif()

execute_process(
    COMMAND ${clean_hook_env_command} "${BENCH}" peeltiming
    --N 16 --bb 2 --target-profile dispatch-v1 --seed-policy raw
    --construction-seed 7 --loss 0.1 --loss-seed 9
    --schedule adversarial --candidate-pmf 0,1
    --candidate-scale identity --warmup-replicates 0 --replicates 4
    --inner-reps 1 --max-overhead 8 --cache-state warm
    --evict-bytes 4096 --context-sha256
    0000000000000000000000000000000000000000000000000000000000000000
    --required-margin 0.05
    RESULT_VARIABLE peeltiming_hard_result
    OUTPUT_VARIABLE peeltiming_hard_out
    ERROR_VARIABLE peeltiming_hard_err
    TIMEOUT 10)
if(NOT peeltiming_hard_result EQUAL 0 OR peeltiming_hard_err OR
   NOT peeltiming_hard_out MATCHES
       "schedule=adversarial,loss_model=packet-schedule-v1,trace_encoding=wirehair-wh2-peeltiming-loss-trace-v1" OR
   NOT peeltiming_hard_out MATCHES
       "\n# peeltiming_done,complete=1,rows=64,")
    message(FATAL_ERROR
        "hard-schedule peeltiming failed in preferred-suppressed build\n"
        "rc=${peeltiming_hard_result}\nstdout=${peeltiming_hard_out}\n"
        "stderr=${peeltiming_hard_err}")
endif()
endif()

execute_process(
    COMMAND ${clean_hook_env_command} "${BENCH}" preferredtiming
    RESULT_VARIABLE timing_result
    OUTPUT_VARIABLE timing_out
    ERROR_VARIABLE timing_err
    TIMEOUT 10)
if(NOT timing_result EQUAL 1 OR
   NOT timing_err STREQUAL "unknown mode: preferredtiming\n")
    message(FATAL_ERROR
        "preferredtiming leaked into the preferred-suppressed CLI\n"
        "rc=${timing_result}\nstdout=${timing_out}\n"
        "stderr=${timing_err}")
endif()

execute_process(
    COMMAND ${clean_hook_env_command} "${BENCH}" preferredattempt
    RESULT_VARIABLE preferred_result
    OUTPUT_VARIABLE preferred_out
    ERROR_VARIABLE preferred_err
    TIMEOUT 10)
if(NOT preferred_result EQUAL 1 OR
   NOT preferred_err STREQUAL "unknown mode: preferredattempt\n")
    message(FATAL_ERROR
        "preferredattempt did not retain the old unknown-mode behavior\n"
        "rc=${preferred_result}\nstdout=${preferred_out}\n"
        "stderr=${preferred_err}")
endif()

execute_process(
    COMMAND ${clean_hook_env_command} "${BENCH}" selftest
    RESULT_VARIABLE selftest_result
    OUTPUT_VARIABLE selftest_out
    ERROR_VARIABLE selftest_err
    TIMEOUT 10)
# The parser oracle deliberately feeds malformed/retired config spellings
# through the production diagnostic path.  Those expected diagnostics moved
# to stderr when parsing became fail-closed; they are not a self-test failure.
if(NOT selftest_result EQUAL 0 OR
   NOT selftest_out MATCHES "loss boundary oracle: PASS" OR
   NOT selftest_err MATCHES "retired S,H,N1,D2 form")
    message(FATAL_ERROR
        "existing preferred-suppressed CLI command changed\n"
        "rc=${selftest_result}\nstdout=${selftest_out}\n"
        "stderr=${selftest_err}")
endif()
