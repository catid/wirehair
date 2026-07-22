if(NOT DEFINED BENCH)
    message(FATAL_ERROR "BENCH is required")
endif()

execute_process(
    COMMAND "${BENCH}"
    RESULT_VARIABLE usage_result
    OUTPUT_VARIABLE usage_out
    ERROR_VARIABLE usage_err
    TIMEOUT 10)
if(NOT usage_result EQUAL 1 OR
   NOT usage_err MATCHES "^usage: wirehair_v2_bench " OR
   usage_err MATCHES "preferredattempt|preferredtiming")
    message(FATAL_ERROR
        "preferredattempt leaked into no-hooks usage\n"
        "rc=${usage_result}\nstdout=${usage_out}\nstderr=${usage_err}")
endif()

execute_process(
    COMMAND "${BENCH}" preferredtiming
    RESULT_VARIABLE timing_result
    OUTPUT_VARIABLE timing_out
    ERROR_VARIABLE timing_err
    TIMEOUT 10)
if(NOT timing_result EQUAL 1 OR
   NOT timing_err STREQUAL "unknown mode: preferredtiming\n")
    message(FATAL_ERROR
        "preferredtiming leaked into the production/no-hooks CLI\n"
        "rc=${timing_result}\nstdout=${timing_out}\n"
        "stderr=${timing_err}")
endif()

execute_process(
    COMMAND "${BENCH}" preferredattempt
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
    COMMAND "${BENCH}" selftest
    RESULT_VARIABLE selftest_result
    OUTPUT_VARIABLE selftest_out
    ERROR_VARIABLE selftest_err
    TIMEOUT 10)
if(NOT selftest_result EQUAL 0 OR selftest_err OR
   NOT selftest_out MATCHES "loss boundary oracle: PASS")
    message(FATAL_ERROR
        "existing no-hooks CLI command changed\n"
        "rc=${selftest_result}\nstdout=${selftest_out}\n"
        "stderr=${selftest_err}")
endif()
