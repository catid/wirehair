if(NOT DEFINED FINGERPRINT OR FINGERPRINT STREQUAL "")
    message(FATAL_ERROR "FINGERPRINT executable is required")
endif()

# Keep this useful as a standalone script as well as through CTest.  A
# developer shell may carry experiment hooks, but argument-parser tests must
# always exercise the frozen equation configuration.
set(clean_hook_env_command
    "${CMAKE_COMMAND}" -E env
    --unset=WIREHAIR_V2_PEEL_DEGREES
    --unset=WIREHAIR_V2_STAIRCASE_DEGREES
    --unset=WIREHAIR_V2_STAIRCASE_ROW_DEGREES
    --unset=WIREHAIR_V2_STAIRCASE_DEGREE_SCALE
    --unset=WIREHAIR_V2_BAND_TRACKING_X)

function(expect_invalid_max_k value)
    execute_process(
        COMMAND ${clean_hook_env_command}
            "${FINGERPRINT}" --print-goldens
            --max-k "${value}" --profile dispatch-v1
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err
        TIMEOUT 30)
    if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 2 OR
       NOT out STREQUAL "" OR
       NOT err STREQUAL "invalid --max-k value\n")
        message(FATAL_ERROR
            "fingerprint accepted malformed --max-k '${value}'\n"
            "result=${result}\nstdout=${out}\nstderr=${err}")
    endif()
endfunction()

foreach(value IN ITEMS "2junk" "2.0" "+2" "-2" " 2" "2 ")
    expect_invalid_max_k("${value}")
endforeach()

execute_process(
    COMMAND ${clean_hook_env_command}
        "${FINGERPRINT}" --print-goldens
        --max-k 2 --profile dispatch-v1
    RESULT_VARIABLE result
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err
    TIMEOUT 30)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0 OR
   NOT out MATCHES
       "static const char kDispatchV1AllKFingerprint\\[\\] =[\r\n ]+\"[0-9a-f]+\";" OR
   NOT err MATCHES "computing dispatch-v1 fingerprint \\(K=2\\.\\.2\\)")
    message(FATAL_ERROR
        "fingerprint rejected valid --max-k 2\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
