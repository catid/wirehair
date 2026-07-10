if(NOT DEFINED BINARY OR NOT EXISTS "${BINARY}")
    message(FATAL_ERROR "missing gen_tables binary: ${BINARY}")
endif()

foreach(run IN ITEMS first second)
    execute_process(
        COMMAND "${BINARY}" --no-benchmarks --heavy-trials 4096
        RESULT_VARIABLE result
        OUTPUT_VARIABLE ${run}_output
        ERROR_VARIABLE ${run}_error)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR
            "gen_tables ${run} run failed (${result})\n${${run}_error}")
    endif()
endforeach()

if(NOT first_output STREQUAL second_output)
    message(FATAL_ERROR "gen_tables heavy output is not deterministic")
endif()

foreach(required IN ITEMS
        "Shipped Cauchy matrix reproduced from seed = 2318331135281"
        "* Empirical perturbation singular rate: 14 / 4096 (expected ~1/256 = 16)"
        "* Informational measurement only; not a reliability guarantee")
    string(FIND "${first_output}" "${required}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR "missing heavy generator evidence: ${required}")
    endif()
endforeach()

string(FIND "${first_output}" "Tests failed" failure_position)
if(NOT failure_position EQUAL -1)
    message(FATAL_ERROR "gen_tables reported a failed regression gate")
endif()
