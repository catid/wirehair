if(NOT DEFINED BENCH)
    message(FATAL_ERROR "BENCH is required")
endif()

function(run_bench result_var out_var err_var)
    execute_process(
        COMMAND "${BENCH}" ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err
        TIMEOUT 90)
    set(${result_var} "${result}" PARENT_SCOPE)
    set(${out_var} "${out}" PARENT_SCOPE)
    set(${err_var} "${err}" PARENT_SCOPE)
endfunction()

function(reject_sanitizer output context)
    if("${output}" MATCHES
       "AddressSanitizer|LeakSanitizer|UndefinedBehaviorSanitizer|MemorySanitizer|ThreadSanitizer|runtime error:")
        message(FATAL_ERROR
            "sanitizer failure during ${context}\n${output}")
    endif()
endfunction()

function(expect_failure pattern)
    run_bench(result out err ${ARGN})
    if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1)
        message(FATAL_ERROR
            "expected exit 1, got '${result}': ${ARGN}\n"
            "stdout=${out}\nstderr=${err}")
    endif()
    if(NOT err MATCHES "${pattern}")
        message(FATAL_ERROR
            "missing diagnostic '${pattern}': ${ARGN}\nstderr=${err}")
    endif()
    reject_sanitizer("${out}${err}" "expected failure: ${ARGN}")
endfunction()

function(expect_success pattern)
    run_bench(result out err ${ARGN})
    if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0)
        message(FATAL_ERROR
            "expected exit 0, got '${result}': ${ARGN}\n"
            "stdout=${out}\nstderr=${err}")
    endif()
    if(NOT out MATCHES "${pattern}")
        message(FATAL_ERROR
            "missing output '${pattern}': ${ARGN}\nstdout=${out}")
    endif()
    reject_sanitizer("${out}${err}" "expected success: ${ARGN}")
endfunction()

# Trial count boundaries, including the old uint16 narrowing boundary.
expect_failure("trials must be" seedtable --N 2 --bb-list 1
    --peel-candidates 1 --trials 0)
foreach(trials IN ITEMS 1 65535 65536 1000000)
    expect_success(",${trials},${trials}," seedtable --N 2 --bb-list 1
        --peel-candidates 1 --trials ${trials})
endforeach()
expect_failure("trials must be" seedtable --N 2 --bb-list 1
    --peel-candidates 1 --trials 1000001)
expect_failure("bad --trials value" seedtable --N 2 --bb-list 1
    --peel-candidates 1 --trials 4294967296)

# Every loss-accepting mode rejects values above the shared boundary.
expect_failure("loss must be" compare --nlo 2 --nhi 2 --trials 1
    --bb-list 1 --max-message-mib 1 --loss 0.9900001)
expect_failure("loss must be" densecheck --N 2 --bb 1 --candidates 1
    --trials 1 --loss 0.9900001)
expect_failure("loss must be" densetune --N 2 --bb-list 1 --candidates 1
    --trials 1 --loss 0.9900001)
expect_failure("loss must be" densecount --N 2 --bb-list 1 --deltas 0
    --trials 1 --loss 0.9900001)
expect_failure("loss must be" densegrid --N 2 --bb-list 1 --deltas 0
    --candidates 1 --trials 1 --loss 0.9900001)

# Parser and interval boundaries.
foreach(loss IN ITEMS -0.1 1.0 nan inf malformed)
    expect_failure("loss" densecheck --N 2 --bb 1 --candidates 1
        --trials 1 --loss ${loss})
    expect_failure("loss" compare --nlo 2 --nhi 2 --trials 1
        --bb-list 1 --max-message-mib 1 --loss ${loss})
    expect_failure("loss" densetune --N 2 --bb-list 1 --candidates 1
        --trials 1 --loss ${loss})
    expect_failure("loss" densecount --N 2 --bb-list 1 --deltas 0
        --trials 1 --loss ${loss})
    expect_failure("loss" densegrid --N 2 --bb-list 1 --deltas 0
        --candidates 1 --trials 1 --loss ${loss})
endforeach()
foreach(loss IN ITEMS 0 0.37)
    expect_success("loss=${loss}" densecheck --N 2 --bb 1 --candidates 1
        --trials 1 --loss ${loss})
    expect_success("loss=${loss}" compare --nlo 2 --nhi 2 --trials 1
        --bb-list 1 --max-message-mib 1 --loss ${loss})
    expect_success("loss=${loss}" densetune --N 2 --bb-list 1
        --candidates 1 --trials 1 --loss ${loss})
    expect_success("loss=${loss}" densecount --N 2 --bb-list 1
        --deltas 0 --trials 1 --loss ${loss})
    expect_success("loss=${loss}" densegrid --N 2 --bb-list 1
        --deltas 0 --candidates 1 --trials 1 --loss ${loss})
endforeach()
expect_success("loss=0.98999999999999999" densecheck --N 2 --bb 1
    --candidates 1 --trials 1 --loss 0.99)
expect_success("loss=0.98999999999999999" compare --nlo 2 --nhi 2
    --trials 1 --bb-list 1 --max-message-mib 1 --loss 0.99)
expect_success("loss=0.98999999999999999" densetune --N 2 --bb-list 1
    --candidates 1 --trials 1 --loss 0.99)
expect_success("loss=0.98999999999999999" densecount --N 2 --bb-list 1
    --deltas 0 --trials 1 --loss 0.99)
expect_success("loss=0.98999999999999999" densegrid --N 2 --bb-list 1
    --deltas 0 --candidates 1 --trials 1 --loss 0.99)

# Valid one-trial smoke for all remaining modes.
expect_success("loss boundary oracle: PASS" selftest)
expect_success("# compare:" compare --nlo 2 --nhi 2 --trials 1
    --bb-list 8 --max-message-mib 1 --loss 0)
expect_success("v2_precode" compare --nlo 64 --nhi 64 --trials 1
    --bb-list 8 --max-message-mib 1 --loss 0 --precode)
expect_success("# densetune:" densetune --N 2 --bb-list 8 --candidates 1
    --trials 1 --loss 0)
expect_success("# densecount:" densecount --N 2 --bb-list 8 --deltas 0
    --trials 1 --loss 0)
expect_success("# densegrid:" densegrid --N 2 --bb-list 8 --deltas 0
    --candidates 1 --trials 1 --loss 0)
expect_success("# peelcost:" peelcost --N 2 --bb-list 8 --trials 1
    --structures lt_m1_c16 --precode dense --overhead 0)
expect_success("# precodefail:" precodefail --N 64 --bb-list 8
    --overhead 0,1 --trials 4 --threads 2 --loss 0.1)
expect_failure("thread launch failed" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 4 --threads 2 --loss 0.1
    --fail-thread-launch-after 1)

# Invalid experiment grids must fail before emitting a result header.
run_bench(result out err peelcost --N 2 --bb-list 8 --trials 1
    --structures lt_m1_c16 --precode dense --overhead -1)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1 OR
    NOT out STREQUAL "" OR
    NOT err MATCHES "overhead must be non-negative")
    message(FATAL_ERROR
        "peelcost invalid-grid failure emitted partial output\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "invalid peelcost grid")

# Strict message caps and widened size checks.
expect_failure("cap cannot accommodate" compare --nlo 2 --nhi 2 --trials 1
    --bb-list 524289 --max-message-mib 1 --loss 0)
expect_success("# compare:" compare --nlo 2 --nhi 2 --trials 1
    --bb-list 524288 --max-message-mib 1 --loss 0)
expect_failure("working set" densecheck --N 64000 --bb 2147483647
    --candidates 1 --trials 1 --loss 0)

if(UNIX AND NOT SKIP_CONSTRAINED)
    set(constrained_commands
        "densecheck --N 2 --bb 100000000 --candidates 1 --trials 1 --loss 0"
        "densetune --N 2 --bb-list 100000000 --candidates 1 --trials 1 --loss 0"
        "densecount --N 2 --bb-list 100000000 --deltas 0 --trials 1 --loss 0"
        "densegrid --N 2 --bb-list 100000000 --deltas 0 --candidates 1 --trials 1 --loss 0")
    foreach(command IN LISTS constrained_commands)
        execute_process(
            COMMAND bash -c "ulimit -v 65536; exec '${BENCH}' ${command}"
            RESULT_VARIABLE result
            OUTPUT_VARIABLE out
            ERROR_VARIABLE err
            TIMEOUT 30)
        if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1 OR
            NOT out STREQUAL "" OR
            NOT err MATCHES "working set cannot be allocated")
            message(FATAL_ERROR
                "constrained-memory failure was not clean: ${command}\n"
                "result=${result}\nstdout=${out}\nstderr=${err}")
        endif()
        reject_sanitizer("${out}${err}" "constrained-memory failure: ${command}")
    endforeach()

    execute_process(
        COMMAND bash -c
            "ulimit -v 65536; exec '${BENCH}' peelcost --N 2 --bb-list 8 --trials 4294967295 --structures lt_m1_c16 --precode dense --overhead 0"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err
        TIMEOUT 30)
    if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1 OR
        NOT out STREQUAL "" OR
        NOT err MATCHES "evaluation storage cannot be allocated")
        message(FATAL_ERROR
            "constrained peelcost failure was not clean\n"
            "result=${result}\nstdout=${out}\nstderr=${err}")
    endif()
    reject_sanitizer("${out}${err}" "constrained peelcost failure")

    execute_process(
        COMMAND bash -c
            "ulimit -v 65536; exec '${BENCH}' compare --nlo 2 --nhi 2 --trials 1 --bb-list 100000000 --max-message-mib 64 --loss 0"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err
        TIMEOUT 30)
    if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1 OR
        NOT out STREQUAL "" OR
        NOT err MATCHES "cap cannot accommodate")
        message(FATAL_ERROR
            "constrained compare failure was not clean\n"
            "result=${result}\nstdout=${out}\nstderr=${err}")
    endif()
    reject_sanitizer("${out}${err}" "constrained compare failure")
endif()
