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

# Candidate work accounting preserves the raw request while reporting the
# bounded, duplicate-free byte-domain work actually completed.
expect_success(",300,256,256,1,1," seedtable --N 2 --bb-list 1
    --peel-candidates 300 --trials 1)
expect_success(",0,1,1,1,1," seedtable --N 2 --bb-list 1
    --peel-candidates 0 --trials 1)
expect_success("2,1,300,256,256," densetune --N 2 --bb-list 1
    --candidates 300 --trials 1 --loss 0)
expect_success("2,1,0,1,1," densetune --N 2 --bb-list 1
    --candidates 0 --trials 1 --loss 0)

# Every loss-accepting mode rejects values above the shared boundary.
expect_failure("loss must be" compare --nlo 2 --nhi 2 --trials 1
    --bb-list 1 --max-message-mib 1 --loss 0.9900001)
expect_failure("loss must be" precodecheck --N 2 --bb-list 1 --trials 1
    --loss 0.9900001)
expect_failure("loss must be" precodefail --N 64 --bb-list 1 --overhead 0
    --trials 1 --threads 1 --loss 0.9900001)
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
    expect_failure("loss" precodecheck --N 2 --bb-list 1 --trials 1
        --loss ${loss})
    expect_failure("loss" precodefail --N 64 --bb-list 1 --overhead 0
        --trials 1 --threads 1 --loss ${loss})
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
    expect_success("loss=${loss}" precodecheck --N 2 --bb-list 1
        --trials 1 --loss ${loss})
    expect_success("loss=${loss}" precodefail --N 64 --bb-list 1
        --overhead 0 --trials 1 --threads 1 --loss ${loss})
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
expect_success("loss=0.98999999999999999" precodecheck --N 2 --bb-list 1
    --trials 1 --loss 0.99)
expect_success("loss=0.98999999999999999" precodefail --N 64 --bb-list 1
    --overhead 0 --trials 1 --threads 1 --loss 0.99)
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
expect_success("64,8,1,1,0,0,0,0,0," precodecheck --N 64 --bb-list 8
    --trials 1 --loss 0)
expect_success("v2_precode[ ]+8[ ]+1[ ]+0" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode)
expect_success("precode_profile=certified" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode)
expect_success("precode_profile_handoff=encoder-selected-v1" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode)
expect_success("v2_mixed[ ]+8[ ]+1[ ]+0" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed)
expect_success("v2_mixed_cached[ ]+8[ ]+1[ ]+0" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode-cache --precode-profile mixed)
expect_success("mixed_period=64" compare --nlo 64 --nhi 64 --trials 1
    --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-period 64)
expect_failure("mixed experiment flags require a mixed precode profile" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile certified --mixed-period 64)
expect_failure("--mixed-period must be in" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-period 11)
expect_failure("--mixed-period must be in" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-period 245)
expect_success("mixed_geometry=shared-x" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-period 64 --mixed-geometry shared-x)
expect_success("mixed_residue_skew=14" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-gf16-rows 4 --mixed-period 29
    --mixed-geometry shared-x --mixed-residue-skew 14)
expect_success("mixed_residue_schedule=ramp" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-gf16-rows 4 --mixed-period 28
    --mixed-geometry shared-x --mixed-residue-schedule ramp)
expect_success("mixed_residue_schedule=hashed mixed_residue_hash_seed=0x7"
    compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-gf16-rows 4 --mixed-period 28
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7)
expect_failure("--mixed-residue-skew must be a corner-preserving" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode --precode-profile mixed --mixed-gf16-rows 4
    --mixed-period 29 --mixed-geometry shared-x --mixed-residue-skew 16)
expect_failure("nonconstant --mixed-residue-schedule requires" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode --precode-profile mixed --mixed-gf16-rows 4
    --mixed-period 28 --mixed-geometry shared-x --mixed-residue-skew 1
    --mixed-residue-schedule ramp)
expect_failure("--mixed-residue-hash-seed requires hashed" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode --precode-profile mixed --mixed-gf16-rows 4
    --mixed-period 28 --mixed-geometry shared-x
    --mixed-residue-schedule ramp --mixed-residue-hash-seed 7)
expect_success("mixed_gf16_rows=3.*mixed_geometry=shared-x" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode --precode-profile mixed --mixed-gf16-rows 3
    --mixed-period 64 --mixed-geometry shared-x)
expect_failure("--mixed-gf16-rows must be in" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-gf16-rows 1)
expect_failure("--mixed-period must be in" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-gf16-rows 3 --mixed-period 12)
expect_failure("unknown --mixed-geometry" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed --mixed-geometry unknown)
expect_failure("mixed experiment flags require a mixed precode profile" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0 --precode --precode-profile certified
    --mixed-geometry shared-x)
expect_success("v2_precode[ ]+17[ ]+1[ ]+0" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 17 --max-message-mib 1 --loss 0 --precode
    --precode-profile certified)
expect_success("v2_cached[ ]+8[ ]+1[ ]+0" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode-cache)
expect_success("encoder_cache=0 decoder_cache=1" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode-decoder-cache)
run_bench(result out err compare --nlo 64 --nhi 64 --trials 2
    --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile both --trial-details)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 0 OR
    NOT out MATCHES "precode_profile=both" OR
    NOT out MATCHES "v2_precode[ ]+8[ ]+2[ ]+0" OR
    NOT out MATCHES "v2_mixed[ ]+8[ ]+2[ ]+0" OR
    NOT out MATCHES
        "paired_trial:.*precode_profile=certified.*precode_ok=1" OR
    NOT out MATCHES "paired_trial:.*precode_profile=mixed.*precode_ok=1")
    message(FATAL_ERROR
        "paired certified/mixed precode comparison failed\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "paired certified/mixed precode comparison")
expect_failure("unknown --precode-profile" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0 --precode
    --precode-profile unknown)
expect_failure("--precode-profile requires" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0
    --precode-profile mixed)
run_bench(result out err compare --nlo 64 --nhi 64 --trials 1
    --bb-list 17 --max-message-mib 1 --loss 0 --precode
    --precode-profile mixed)
if(NOT result MATCHES "^-?[0-9]+$" OR NOT result EQUAL 1 OR
    NOT out STREQUAL "" OR
    NOT err MATCHES "mixed precode profile requires even block bytes")
    message(FATAL_ERROR
        "mixed odd-byte rejection emitted partial output or wrong status\n"
        "result=${result}\nstdout=${out}\nstderr=${err}")
endif()
reject_sanitizer("${out}${err}" "mixed odd-byte rejection")
expect_failure("mixed precode profile requires even block bytes" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8,17 --max-message-mib 1
    --loss 0 --precode --precode-profile both)
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
string(CONCAT precodefail_paired_pattern
    "precodefail_paired:.*completion=mixed.*mix2_fail=.*mix3_fail=.*"
    "both_fail=.*mix2_only=.*mix3_only=.*both_success=.*mcnemar_p=.*"
    "mix2_seed_attempt=.*mix3_seed_attempt=.*mix2_wilson95=.*"
    "mix3_wilson95=")
expect_success("${precodefail_paired_pattern}" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 4 --threads 2 --loss 0.1
    --mix-count 2,3 --completion mixed --payload-e2e)
expect_success("mixed_period=64 mixed_gf16_rows=2 mixed_geometry=shared-x.*full_payload_solve=1"
    precodefail
    --N 64 --bb-list 1280 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-period 64 --mixed-geometry shared-x
    --full-payload-solve)
expect_success("mixed_period=64 mixed_gf16_rows=3 mixed_geometry=shared-x"
    precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mixed-gf16-rows 3 --mixed-period 64
    --mixed-geometry shared-x)
expect_success("mixed_period=32 mixed_gf16_rows=4 mixed_geometry=shared-x"
    precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 32
    --mixed-geometry shared-x)
expect_success("mixed_residue_skew=14" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 29
    --mixed-geometry shared-x --mixed-residue-skew 14)
expect_success("mixed_residue_schedule=ramp" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 28
    --mixed-geometry shared-x --mixed-residue-schedule ramp)
expect_success("mixed_residue_schedule=hashed mixed_residue_hash_seed=0x7"
    precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 2 --threads 2 --loss 0.1
    --completion mixed --mixed-gf16-rows 4 --mixed-period 28
    --mixed-geometry shared-x --mixed-residue-schedule hashed
    --mixed-residue-hash-seed 7)
expect_failure("--mixed-gf16-rows must be in" precodefail --N 64
    --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --mixed-gf16-rows 5)
expect_failure("--mixed-period must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1 --completion mixed
    --mixed-gf16-rows 3 --mixed-period 12)
expect_failure("mixed experiment flags require --completion mixed" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --mixed-period 64)
expect_failure("--mixed-period must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1 --completion mixed
    --mixed-period 11)
expect_failure("--mixed-period must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1 --completion mixed
    --mixed-period 245)
expect_failure("unknown --mixed-geometry" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1 --completion mixed
    --mixed-geometry unknown)
expect_failure("mixed experiment flags require --completion mixed" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --mixed-geometry shared-x)
expect_failure("unknown --completion" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion unknown)
expect_failure("mixed completion requires even block bytes" precodefail --N 64
    --bb-list 7 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed)
expect_failure("mixed completion requires periodic heavy family" precodefail
    --N 64 --bb-list 8 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --completion mixed --heavy-family hashed)
expect_failure("working set exceeds|message dimensions are unsupported" precodefail --N 64000
    --bb-list 2147483647 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --payload-e2e)
expect_failure("working set exceeds|message dimensions are unsupported" precodefail --N 64000
    --bb-list 2200000 --overhead 0 --trials 1 --threads 1 --loss 0.1
    --payload-e2e)

# Common packet schedules are accepted by both E2E comparison entry points.
foreach(schedule IN ITEMS iid burst permutation systematic-first repair-only
        adversarial)
    expect_success("schedule=${schedule}" compare --nlo 64 --nhi 64
        --trials 1 --bb-list 8 --max-message-mib 1 --loss 0.1
        --schedule ${schedule} --precode)
    expect_success("schedule=${schedule}" precodecheck --N 64 --bb-list 8
        --trials 1 --loss 0.1 --schedule ${schedule})
endforeach()
expect_success("paired_trial: schedule=burst" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --loss 0.1
    --schedule burst --precode --trial-details)
expect_success("cached_ok=1.*cached_oh=[0-9]+.*cached_delta=" compare
    --nlo 64 --nhi 64 --trials 1 --bb-list 8 --max-message-mib 1
    --loss 0.1 --schedule iid --precode-cache --trial-details)
expect_success("precode_trial: schedule=permutation" precodecheck --N 64
    --bb-list 8 --trials 1 --loss 0.1 --schedule permutation
    --trial-details)
# High burst loss used to clamp its start probability to one and drop forever.
# Exercise both consumers at the accepted loss boundary.
expect_success("schedule=burst" compare --nlo 8 --nhi 8 --trials 1
    --bb-list 8 --max-message-mib 1 --loss 0.99 --schedule burst --precode)
expect_success("schedule=burst" precodecheck --N 64 --bb-list 8
    --trials 1 --loss 0.99 --schedule burst)
expect_failure("unknown compare schedule" compare --nlo 64 --nhi 64
    --trials 1 --bb-list 8 --max-message-mib 1 --schedule unknown)
expect_failure("unknown precodecheck schedule" precodecheck --N 64
    --bb-list 8 --trials 1 --schedule unknown)
expect_success("hashed" precodefail --N 64 --bb-list 1
    --overhead 0 --trials 1 --threads 1 --loss 0.1
    --heavy-family periodic,hashed)
expect_failure("unknown --heavy-family" precodefail --N 64 --bb-list 1
    --overhead 0 --trials 1 --threads 1 --loss 0.1
    --heavy-family unknown)

# Thread parsing and the partially-launched worker cleanup path.
expect_failure("bad --threads value" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads malformed --loss 0.1)
expect_failure("threads must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 0 --loss 0.1)
expect_failure("threads must be in" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 257 --loss 0.1)
expect_failure("bad --threads value" precodefail --N 64 --bb-list 8
    --overhead 0 --trials 1 --threads 4294967296 --loss 0.1)
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
