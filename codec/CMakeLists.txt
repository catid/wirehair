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

    add_executable(wirehair_v2_precode_test PrecodeTest.cpp)
    wirehair_apply_project_options(wirehair_v2_precode_test)
    target_link_libraries(wirehair_v2_precode_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_precode_test COMMAND wirehair_v2_precode_test)

    add_executable(wirehair_v2_precode_encode_test PrecodeEncodeTest.cpp)
    wirehair_apply_project_options(wirehair_v2_precode_encode_test)
    target_link_libraries(wirehair_v2_precode_encode_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_precode_encode_test COMMAND wirehair_v2_precode_encode_test)

    add_executable(wirehair_v2_precode_solve_test PrecodeSolveTest.cpp)
    wirehair_apply_project_options(wirehair_v2_precode_solve_test)
    target_link_libraries(wirehair_v2_precode_solve_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_precode_solve_test COMMAND wirehair_v2_precode_solve_test)

    add_executable(wirehair_v2_precode_roundtrip_test PrecodeRoundTripTest.cpp)
    wirehair_apply_project_options(wirehair_v2_precode_roundtrip_test)
    target_link_libraries(wirehair_v2_precode_roundtrip_test wirehair_v2_policy)
    add_test(NAME wirehair_v2_precode_roundtrip_test COMMAND wirehair_v2_precode_roundtrip_test)

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
if(BUILD_TESTS AND WIREHAIR_BUILD_BENCHMARKS)
    add_test(
        NAME wirehair_v2_bench_cli_test
        COMMAND ${CMAKE_COMMAND}
            -DBENCH=$<TARGET_FILE:wirehair_v2_bench>
            -DSKIP_CONSTRAINED=${_wirehair_v2_sanitized}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/V2BenchCliTest.cmake)
    set_tests_properties(wirehair_v2_bench_cli_test PROPERTIES TIMEOUT 120)
endif()
