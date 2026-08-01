add_library(wirehair_v2_policy STATIC
    WirehairV2Codec.cpp
    WirehairV2Codec.h
    WirehairV2Peel.cpp
    WirehairV2Peel.h
    WirehairV2Plan.cpp
    WirehairV2Plan.h
    WirehairV2Policy.cpp
    WirehairV2Policy.h
    WirehairV2Precode.cpp
    WirehairV2Precode.h
    WirehairV2PrecodeDecode.cpp
    WirehairV2PrecodeDecode.h
    WirehairV2PrecodeEncode.cpp
    WirehairV2PrecodeEncode.h
    WirehairV2Seeds.cpp
    WirehairV2Seeds.h
    WirehairV2Solve.cpp
    WirehairV2Solve.h
)

target_include_directories(wirehair_v2_policy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(wirehair_v2_policy PUBLIC ${WIREHAIR_INTERNAL_LIBRARY})
wirehair_apply_project_options(wirehair_v2_policy)

if (BUILD_TESTS)
    target_compile_definitions(wirehair_v2_policy
        PUBLIC WIREHAIR_V2_ENABLE_TEST_HOOKS=1)

    set(_wirehair_v2_compile_flags
        "${CMAKE_CXX_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
    if(CMAKE_CONFIGURATION_TYPES)
        foreach(_wirehair_v2_config IN LISTS CMAKE_CONFIGURATION_TYPES)
            string(TOUPPER "${_wirehair_v2_config}"
                _wirehair_v2_config_upper)
            string(APPEND _wirehair_v2_compile_flags
                " ${CMAKE_CXX_FLAGS_${_wirehair_v2_config_upper}}"
                " ${CMAKE_EXE_LINKER_FLAGS_${_wirehair_v2_config_upper}}")
        endforeach()
    else()
        string(TOUPPER "${CMAKE_BUILD_TYPE}" _wirehair_v2_build_type)
        string(APPEND _wirehair_v2_compile_flags
            " ${CMAKE_CXX_FLAGS_${_wirehair_v2_build_type}}"
            " ${CMAKE_EXE_LINKER_FLAGS_${_wirehair_v2_build_type}}")
    endif()
    set(_wirehair_v2_sanitized FALSE)
    if(_wirehair_v2_compile_flags MATCHES "(^|[ ;])-fsanitize=")
        set(_wirehair_v2_sanitized TRUE)
    endif()

    add_executable(wirehair_v2_policy_test PolicyTest.cpp)
    wirehair_apply_project_options(wirehair_v2_policy_test)
    target_link_libraries(wirehair_v2_policy_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_policy_test COMMAND wirehair_v2_policy_test)

    add_executable(wirehair_v2_noinit_test V2NoInitTest.cpp)
    wirehair_apply_project_options(wirehair_v2_noinit_test)
    target_link_libraries(wirehair_v2_noinit_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_noinit_test COMMAND wirehair_v2_noinit_test)

    add_executable(wirehair_v2_transactional_test V2TransactionalTest.cpp)
    wirehair_apply_project_options(wirehair_v2_transactional_test)
    target_link_libraries(wirehair_v2_transactional_test
        wirehair_v2_policy Threads::Threads)
    add_test(NAME wirehair_v2_transactional_test
        COMMAND wirehair_v2_transactional_test)
    set_tests_properties(wirehair_v2_transactional_test
        PROPERTIES LABELS "state-machine;transactional;v2" TIMEOUT 180)

    add_executable(wirehair_v2_precode_test PrecodeTest.cpp)
    wirehair_apply_project_options(wirehair_v2_precode_test)
    target_link_libraries(wirehair_v2_precode_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_precode_test COMMAND wirehair_v2_precode_test)

    add_executable(wirehair_v2_precode_encode_test PrecodeEncodeTest.cpp)
    wirehair_apply_project_options(wirehair_v2_precode_encode_test)
    target_link_libraries(wirehair_v2_precode_encode_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_precode_encode_test COMMAND wirehair_v2_precode_encode_test)

    add_executable(wirehair_v2_precode_decode_test PrecodeDecodeTest.cpp)
    wirehair_apply_project_options(wirehair_v2_precode_decode_test)
    target_link_libraries(wirehair_v2_precode_decode_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_precode_decode_test COMMAND wirehair_v2_precode_decode_test)

    add_executable(wirehair_v2_precode_solve_test
        PrecodeSolveTest.cpp
        V2TinyDenseOracle.cpp
        V2TinyDenseOracle.h)
    wirehair_apply_project_options(wirehair_v2_precode_solve_test)
    target_link_libraries(wirehair_v2_precode_solve_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_precode_solve_test
        COMMAND wirehair_v2_precode_solve_test)

    add_executable(wirehair_v2_precode_roundtrip_test PrecodeRoundTripTest.cpp)
    wirehair_apply_project_options(wirehair_v2_precode_roundtrip_test)
    target_link_libraries(wirehair_v2_precode_roundtrip_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_precode_roundtrip_test
        COMMAND wirehair_v2_precode_roundtrip_test)

    add_executable(wirehair_v2_precode_seed_selection_test
        PrecodeSeedSelectionTest.cpp)
    wirehair_apply_project_options(wirehair_v2_precode_seed_selection_test)
    target_link_libraries(wirehair_v2_precode_seed_selection_test
        wirehair_v2_policy)
    if(_wirehair_v2_sanitized)
        set(_wirehair_v2_seed_selection_args --sanitizer)
    else()
        set(_wirehair_v2_seed_selection_args "")
    endif()
    add_test(NAME wirehair_v2_precode_seed_selection_test
        COMMAND wirehair_v2_precode_seed_selection_test
            ${_wirehair_v2_seed_selection_args})
    set_tests_properties(wirehair_v2_precode_seed_selection_test
        PROPERTIES TIMEOUT 180)

    add_executable(wirehair_v2_profile_test V2ProfileTest.cpp)
    wirehair_apply_project_options(wirehair_v2_profile_test)
    target_link_libraries(wirehair_v2_profile_test wirehair)
    if(WIN32 AND WIREHAIR_PRIMARY_VARIANT STREQUAL "SHARED")
        # Windows has no build-tree RPATH.  This executable lives in the
        # codec configuration directory, while the primary DLL lives in the
        # root configuration directory, so stage the runtime beside it.
        add_custom_command(TARGET wirehair_v2_profile_test POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:wirehair>
                $<TARGET_FILE_DIR:wirehair_v2_profile_test>)
    endif()
    add_test(NAME wirehair_v2_profile_test COMMAND wirehair_v2_profile_test)
    set_tests_properties(wirehair_v2_profile_test PROPERTIES
        LABELS "v2;profile;serialization;golden;e2e;abi"
        TIMEOUT 120)

    set(_wirehair_v2_fuzz_common
        fuzz/V2FuzzDriver.cpp
        fuzz/V2FuzzDriver.h)
    add_executable(wirehair_v2_profile_fuzz
        fuzz/V2ProfileFuzz.cpp
        V2TinyDenseOracle.cpp
        V2TinyDenseOracle.h
        ${_wirehair_v2_fuzz_common})
    target_compile_definitions(wirehair_v2_profile_fuzz PRIVATE
        "WIREHAIR_FUZZ_TARGET_NAME=\"profile\""
        "WIREHAIR_FUZZ_CORPUS_MANIFEST=\"${CMAKE_CURRENT_SOURCE_DIR}/fuzz/corpus/profile/manifest.txt\"")
    target_link_libraries(wirehair_v2_profile_fuzz wirehair_v2_policy)
    wirehair_apply_project_options(wirehair_v2_profile_fuzz)

    add_executable(wirehair_v2_stateful_fuzz
        fuzz/V2StatefulFuzz.cpp
        ${_wirehair_v2_fuzz_common})
    target_compile_definitions(wirehair_v2_stateful_fuzz PRIVATE
        "WIREHAIR_FUZZ_TARGET_NAME=\"stateful\""
        "WIREHAIR_FUZZ_CORPUS_MANIFEST=\"${CMAKE_CURRENT_SOURCE_DIR}/fuzz/corpus/stateful/manifest.txt\"")
    target_link_libraries(wirehair_v2_stateful_fuzz wirehair_v2_policy)
    wirehair_apply_project_options(wirehair_v2_stateful_fuzz)

    if(WIREHAIR_ENABLE_LIBFUZZER)
        if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            message(FATAL_ERROR
                "WIREHAIR_ENABLE_LIBFUZZER requires Clang")
        endif()
        foreach(_wirehair_v2_fuzz_target IN ITEMS
                wirehair_v2_profile_fuzz wirehair_v2_stateful_fuzz)
            target_compile_definitions(${_wirehair_v2_fuzz_target}
                PRIVATE WIREHAIR_ENABLE_LIBFUZZER=1)
            target_compile_options(${_wirehair_v2_fuzz_target}
                PRIVATE -fsanitize=fuzzer)
            target_link_options(${_wirehair_v2_fuzz_target}
                PRIVATE -fsanitize=fuzzer)
        endforeach()
    else()
        set(_wirehair_v2_fuzz_artifacts
            "${CMAKE_BINARY_DIR}/Testing/Temporary/v2-fuzz-artifacts")
        file(MAKE_DIRECTORY "${_wirehair_v2_fuzz_artifacts}")
        add_test(NAME wirehair_v2_profile_fuzz_test
            COMMAND wirehair_v2_profile_fuzz
                --mutations 10000
                --max-seconds 55
                --seed 0x243f6a8885a308d3
                --trace-current
                --artifact-dir "${_wirehair_v2_fuzz_artifacts}")
        add_test(NAME wirehair_v2_stateful_fuzz_test
            COMMAND wirehair_v2_stateful_fuzz
                --mutations 10000
                --max-seconds 55
                --seed 0x13198a2e03707344
                --trace-current
                --artifact-dir "${_wirehair_v2_fuzz_artifacts}")
        set_tests_properties(
            wirehair_v2_profile_fuzz_test
            wirehair_v2_stateful_fuzz_test
            PROPERTIES
                LABELS "fuzz;sanitizer;v2"
                TIMEOUT 60
                PROCESSORS 1)
    endif()

    if(WIREHAIR_ENABLE_SCHEDULED_TESTS)
        add_test(
            NAME wirehair_v2_precode_decode_large
            COMMAND wirehair_v2_precode_decode_test --large)
        add_test(
            NAME wirehair_v2_precode_decode_large_recovery
            COMMAND wirehair_v2_precode_decode_test --large-recovery)
        set_tests_properties(
            wirehair_v2_precode_decode_large
            wirehair_v2_precode_decode_large_recovery
            PROPERTIES LABELS "scheduled;v2" TIMEOUT 600 RUN_SERIAL TRUE)
    endif()

    if(UNIX AND NOT _wirehair_v2_sanitized)
        add_test(
            NAME wirehair_v2_precode_oom_test
            COMMAND bash -c
                "ulimit -v 65536; exec '$<TARGET_FILE:wirehair_v2_precode_encode_test>' --oom-probe")
        set_tests_properties(wirehair_v2_precode_oom_test PROPERTIES TIMEOUT 30)
    endif()
endif()

if(WIREHAIR_BUILD_BENCHMARKS)
    set(WIREHAIR_BENCH_EXCLUDE "")
else()
    set(WIREHAIR_BENCH_EXCLUDE EXCLUDE_FROM_ALL)
endif()
add_executable(wirehair_v2_bench ${WIREHAIR_BENCH_EXCLUDE}
    WirehairV2Bench.cpp)
target_link_libraries(wirehair_v2_bench wirehair_v2_policy Threads::Threads)
wirehair_apply_project_options(wirehair_v2_bench)

add_executable(wirehair_v2_packet_eval_bench ${WIREHAIR_BENCH_EXCLUDE}
    ${PROJECT_SOURCE_DIR}/experiments/v2_packet_eval_bench.cpp)
target_link_libraries(wirehair_v2_packet_eval_bench wirehair_v2_policy)
wirehair_apply_project_options(wirehair_v2_packet_eval_bench)

add_executable(wirehair_v2_resume_bench ${WIREHAIR_BENCH_EXCLUDE}
    V2ResumeBench.cpp)
target_link_libraries(wirehair_v2_resume_bench wirehair_v2_policy)
wirehair_apply_project_options(wirehair_v2_resume_bench)

if(BUILD_TESTS AND WIREHAIR_BUILD_BENCHMARKS)
    add_test(
        NAME wirehair_v2_bench_cli_test
        COMMAND ${CMAKE_COMMAND}
            -DBENCH=$<TARGET_FILE:wirehair_v2_bench>
            -DSKIP_CONSTRAINED=${_wirehair_v2_sanitized}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/V2BenchCliTest.cmake)
    set_tests_properties(wirehair_v2_bench_cli_test PROPERTIES TIMEOUT 120)
endif()
