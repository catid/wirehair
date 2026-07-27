if(NOT DEFINED BENCH)
    message(FATAL_ERROR "BENCH is required")
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
   usage_err MATCHES "preferredattempt|preferredtiming")
    message(FATAL_ERROR
        "preferredattempt leaked into preferred-suppressed usage\n"
        "rc=${usage_result}\nstdout=${usage_out}\nstderr=${usage_err}")
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
