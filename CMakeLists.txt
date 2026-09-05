cmake_minimum_required(VERSION 3.15)
project(wirehair VERSION 2.0.0 LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckCXXCompilerFlag)

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
endif()
if(NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
find_package(Threads REQUIRED)

# Embed one configure-time source identity in benchmark evidence producers.
# Controllers separately require a clean, stable worktree at execution time.
set(WIREHAIR_WH2_SOURCE_GIT_COMMIT
    "0000000000000000000000000000000000000000")
find_package(Git QUIET)
if(Git_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --verify HEAD^{commit}
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        RESULT_VARIABLE WIREHAIR_WH2_GIT_RESULT
        OUTPUT_VARIABLE WIREHAIR_WH2_GIT_OUTPUT
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(LENGTH "${WIREHAIR_WH2_GIT_OUTPUT}"
        WIREHAIR_WH2_GIT_OUTPUT_LENGTH)
    if(WIREHAIR_WH2_GIT_RESULT EQUAL 0 AND
       WIREHAIR_WH2_GIT_OUTPUT_LENGTH EQUAL 40 AND
       WIREHAIR_WH2_GIT_OUTPUT MATCHES "^[0-9a-f]+$")
        set(WIREHAIR_WH2_SOURCE_GIT_COMMIT
            "${WIREHAIR_WH2_GIT_OUTPUT}")
    endif()
endif()

option(BUILD_SHARED_LIBS "Build Wirehair as a shared library" OFF)
option(BUILD_TESTS "Build and register the correctness tests" ON)
option(BUILD_CODEC_V2
    "Build internal V2 tests, fuzz targets, benchmarks, and tools" ON)
option(MARCH_NATIVE "Tune Wirehair for the build host (non-portable)" OFF)
option(WIREHAIR_BUILD_BOTH "Build both static and shared library variants" OFF)
option(WIREHAIR_STATIC_PIC
    "Build static Wirehair objects as position-independent code" ON)
option(WIREHAIR_BUILD_TOOLS
    "Include offline table and seed generators in the default build" OFF)
option(WIREHAIR_BUILD_BENCHMARKS
    "Include benchmark executables in the default build" OFF)
option(WIREHAIR_ENABLE_SCHEDULED_TESTS
    "Register resource-bounded scheduled large-message profiles" OFF)
option(WIREHAIR_ENABLE_PUBLIC_BORROWED_R1_NATIVE_CONFIG
    "Register the host-specific public borrowed r1 oracle test on CPU 120" OFF)
option(WIREHAIR_ENABLE_LIBFUZZER
    "Build V2 fuzz entry points with Clang libFuzzer main/coverage" OFF)
option(WIREHAIR_ENABLE_TRI_BMQ4_SCHUR_FALSIFIER
    "Build the one-cell tri-bmq4-schur-v0 structural falsifier" OFF)
option(WIREHAIR_STRICT_WARNINGS
    "Treat compiler warnings as errors for Wirehair targets" OFF)

if(WIREHAIR_ENABLE_TRI_BMQ4_SCHUR_FALSIFIER AND
   (NOT BUILD_TESTS OR NOT BUILD_CODEC_V2))
    message(FATAL_ERROR
        "tri-bmq4-schur-v0 falsifier requires BUILD_TESTS and BUILD_CODEC_V2")
endif()

set(WH_LTO "OFF" CACHE STRING "Link-time optimization mode: OFF, ON, AUTO, or THIN")
set_property(CACHE WH_LTO PROPERTY STRINGS OFF ON AUTO THIN)
set(WH_PGO_MODE "OFF" CACHE STRING "Profile-guided optimization mode: OFF, GENERATE, or USE")
set_property(CACHE WH_PGO_MODE PROPERTY STRINGS OFF GENERATE USE)
set(WH_PGO_DIR "${CMAKE_BINARY_DIR}/pgo" CACHE PATH "Profile data directory for WH_PGO_MODE")

set(WIREHAIR_ELF_EXPORT_MAP
    "${PROJECT_SOURCE_DIR}/abi/wirehair.map")
set(WIREHAIR_MACHO_EXPORT_LIST "")
if(APPLE)
    # Keep ELF, PE, and Mach-O on the same versioned C-ABI source of truth.
    # Darwin's linker expects leading underscores in its exported-symbol list.
    file(STRINGS "${WIREHAIR_ELF_EXPORT_MAP}"
        _wirehair_export_lines REGEX "^[ \t]*wirehair_[A-Za-z0-9_]+;")
    if(NOT _wirehair_export_lines)
        message(FATAL_ERROR
            "No public symbols found in ${WIREHAIR_ELF_EXPORT_MAP}")
    endif()
    set(WIREHAIR_MACHO_EXPORT_LIST
        "${CMAKE_CURRENT_BINARY_DIR}/wirehair-macho.exports")
    file(WRITE "${WIREHAIR_MACHO_EXPORT_LIST}" "")
    foreach(_wirehair_export_line IN LISTS _wirehair_export_lines)
        string(STRIP "${_wirehair_export_line}" _wirehair_export_line)
        if(NOT _wirehair_export_line MATCHES
           "^wirehair_[A-Za-z0-9_]+;$")
            message(FATAL_ERROR
                "Malformed public export: '${_wirehair_export_line}'")
        endif()
        string(REGEX REPLACE ";$" "" _wirehair_export_symbol
            "${_wirehair_export_line}")
        file(APPEND "${WIREHAIR_MACHO_EXPORT_LIST}"
            "_${_wirehair_export_symbol}\n")
    endforeach()
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
        "${WIREHAIR_ELF_EXPORT_MAP}")
endif()

set(WIREHAIR_V2_PUBLIC_SOURCE_FILES
        codec/WirehairV2Codec.cpp
        codec/WirehairV2Codec.h
        codec/WirehairV2Peel.cpp
        codec/WirehairV2Peel.h
        codec/WirehairV2Plan.cpp
        codec/WirehairV2Plan.h
        codec/WirehairV2Policy.cpp
        codec/WirehairV2Policy.h
        codec/WirehairV2Precode.cpp
        codec/WirehairV2Precode.h
        codec/WirehairV2PrecodeDecode.cpp
        codec/WirehairV2PrecodeDecode.h
        codec/WirehairV2PrecodeEncode.cpp
        codec/WirehairV2PrecodeEncode.h
        codec/WirehairV2Profile.cpp
        codec/WirehairV2Seeds.cpp
        codec/WirehairV2Seeds.h
        codec/WirehairV2Solve.cpp
        codec/WirehairV2Solve.h)

set(LIB_SOURCE_FILES
        wirehair.cpp
        include/wirehair/wirehair.h
        include/wirehair/wirehair.hpp
        gf256.cpp
        gf256.h
        WirehairCodec.cpp
        WirehairCodec.h
        WirehairCompiler.h
        WirehairEnvironment.h
        WirehairHeavy.h
        WirehairTools.cpp
        WirehairTools.h
        ${WIREHAIR_V2_PUBLIC_SOURCE_FILES}
        )

set(UNIT_TEST_SOURCE_FILES
        test/SiameseTools.cpp
        test/SiameseTools.h
        test/UnitTest.cpp
        )

set(LARGE_MESSAGE_TEST_SOURCE_FILES
        test/LargeMessageTest.cpp
        )

set(GEN_SMALL_DSEEDS
        test/SiameseTools.cpp
        test/SiameseTools.h
        tables/GenerateSmallDenseSeeds.cpp
        )

set(GEN_PEEL_SEEDS
        test/SiameseTools.cpp
        test/SiameseTools.h
        tables/GeneratePeelSeeds.cpp
        )

set(GEN_MOST_DSEEDS
        test/SiameseTools.cpp
        test/SiameseTools.h
        tables/GenerateMostDenseSeeds.cpp
        )

set(GEN_DCOUNTS
        test/SiameseTools.cpp
        test/SiameseTools.h
        tables/GenerateDenseCount.cpp
        )

set(GEN_TABLES
        test/SiameseTools.cpp
        test/SiameseTools.h
        tables/TableGenerator.cpp
        tables/HeavyRowGenerator.cpp
        tables/HeavyRowGenerator.h
        gf256.cpp
        gf256.h
        WirehairCompiler.h
        WirehairEnvironment.h
        )

set(WIREHAIR_NATIVE_COMPILE_OPTION "")
check_cxx_compiler_flag("-march=native" WIREHAIR_HAS_MARCH_NATIVE)
if(MARCH_NATIVE)
    if(NOT WIREHAIR_HAS_MARCH_NATIVE)
        message(FATAL_ERROR
            "MARCH_NATIVE=ON was requested, but ${CMAKE_CXX_COMPILER_ID} "
            "does not support -march=native")
    endif()
    set(WIREHAIR_NATIVE_COMPILE_OPTION "-march=native")
endif()

set(WIREHAIR_OPT_COMPILE_OPTIONS "")
set(WIREHAIR_OPT_LINK_OPTIONS "")
if(NOT MSVC)
    string(TOUPPER "${WH_LTO}" WH_LTO_UC)
    set(WH_LTO_FLAG "")
    if(WH_LTO_UC STREQUAL "ON")
        set(WH_LTO_FLAG "-flto")
    elseif(WH_LTO_UC STREQUAL "AUTO")
        set(WH_LTO_FLAG "-flto=auto")
    elseif(WH_LTO_UC STREQUAL "THIN")
        if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            message(FATAL_ERROR "WH_LTO=THIN requires Clang; ${CMAKE_CXX_COMPILER_ID} does not support -flto=thin")
        endif()
        set(WH_LTO_FLAG "-flto=thin")
    elseif(NOT WH_LTO_UC STREQUAL "OFF")
        message(FATAL_ERROR "Unknown WH_LTO='${WH_LTO}' (expected OFF, ON, AUTO, or THIN)")
    endif()
    if(WH_LTO_FLAG)
        list(APPEND WIREHAIR_OPT_COMPILE_OPTIONS "${WH_LTO_FLAG}")
        list(APPEND WIREHAIR_OPT_LINK_OPTIONS "${WH_LTO_FLAG}")
    endif()

    string(TOUPPER "${WH_PGO_MODE}" WH_PGO_UC)
    if(WH_PGO_UC STREQUAL "GENERATE")
        file(MAKE_DIRECTORY "${WH_PGO_DIR}")
        list(APPEND WIREHAIR_OPT_COMPILE_OPTIONS
            "-fprofile-generate=${WH_PGO_DIR}")
        list(APPEND WIREHAIR_OPT_LINK_OPTIONS
            "-fprofile-generate=${WH_PGO_DIR}")
    elseif(WH_PGO_UC STREQUAL "USE")
        list(APPEND WIREHAIR_OPT_COMPILE_OPTIONS
            "-fprofile-use=${WH_PGO_DIR}"
            -fprofile-correction
            -Wno-missing-profile)
        list(APPEND WIREHAIR_OPT_LINK_OPTIONS
            "-fprofile-use=${WH_PGO_DIR}"
            -fprofile-correction
            -Wno-missing-profile)
    elseif(NOT WH_PGO_UC STREQUAL "OFF")
        message(FATAL_ERROR "Unknown WH_PGO_MODE='${WH_PGO_MODE}' (expected OFF, GENERATE, or USE)")
    endif()
endif()

function(wirehair_apply_project_options target)
    target_compile_features(${target} PRIVATE cxx_std_11)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
        if(WIREHAIR_STRICT_WARNINGS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra)
        if(WIREHAIR_STRICT_WARNINGS)
            target_compile_options(${target} PRIVATE -Wpedantic -Werror)
        endif()
    endif()
    if(WIREHAIR_NATIVE_COMPILE_OPTION)
        target_compile_options(${target} PRIVATE
            "${WIREHAIR_NATIVE_COMPILE_OPTION}")
    endif()
    if(WIREHAIR_OPT_COMPILE_OPTIONS)
        target_compile_options(${target} PRIVATE
            ${WIREHAIR_OPT_COMPILE_OPTIONS})
    endif()
    get_target_property(_wirehair_target_type ${target} TYPE)
    if(WIREHAIR_OPT_LINK_OPTIONS AND
       NOT _wirehair_target_type STREQUAL "STATIC_LIBRARY" AND
       NOT _wirehair_target_type STREQUAL "OBJECT_LIBRARY")
        target_link_options(${target} PRIVATE ${WIREHAIR_OPT_LINK_OPTIONS})
    endif()
endfunction()

function(wirehair_configure_library target variant)
    set_target_properties(${target} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR})
    target_include_directories(${target} PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    # The installed C++ convenience header uses C++11 language/library
    # features.  Publish that requirement so consumers are upgraded even when
    # their project starts from an older CXX_STANDARD.
    target_compile_features(${target} PUBLIC cxx_std_11)
    target_link_libraries(${target} PRIVATE Threads::Threads)
    if(UNIX)
        # The public V2 policy uses libm.  Static archives must publish that
        # dependency to CMake consumers; shared objects resolve it themselves
        # without making downstream applications overlink it.
        if(variant STREQUAL "STATIC")
            target_link_libraries(${target} PUBLIC m)
        else()
            target_link_libraries(${target} PRIVATE m)
        endif()
    endif()
    if(variant STREQUAL "SHARED")
        target_compile_definitions(${target}
            PRIVATE WIREHAIR_BUILDING=1 WIREHAIR_DLL=1
            INTERFACE WIREHAIR_DLL=1)
        if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
            target_link_options(${target} PRIVATE
                "LINKER:--version-script=${WIREHAIR_ELF_EXPORT_MAP}")
            set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS
                "${WIREHAIR_ELF_EXPORT_MAP}")
        elseif(APPLE)
            target_link_options(${target} PRIVATE
                "LINKER:-exported_symbols_list,${WIREHAIR_MACHO_EXPORT_LIST}")
            set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS
                "${WIREHAIR_ELF_EXPORT_MAP}"
                "${WIREHAIR_MACHO_EXPORT_LIST}")
        endif()
    else()
        target_compile_definitions(${target}
            PRIVATE WIREHAIR_BUILDING=1
            PUBLIC WIREHAIR_STATIC=1)
        set_target_properties(${target} PROPERTIES
            POSITION_INDEPENDENT_CODE ${WIREHAIR_STATIC_PIC})
    endif()
    wirehair_apply_project_options(${target})
endfunction()

if(BUILD_SHARED_LIBS)
    set(WIREHAIR_PRIMARY_VARIANT SHARED)
else()
    set(WIREHAIR_PRIMARY_VARIANT STATIC)
endif()

set(WIREHAIR_INSTALL_TARGETS wirehair)
if(WIREHAIR_BUILD_BOTH AND NOT WIN32)
    # ELF/Mach-O static and shared variants have identical core compile state.
    # PIC objects let both variants share one compilation safely.
    add_library(wirehair_objects OBJECT ${LIB_SOURCE_FILES})
    set_target_properties(wirehair_objects PROPERTIES
        POSITION_INDEPENDENT_CODE ON)
    target_include_directories(wirehair_objects PRIVATE
        ${PROJECT_SOURCE_DIR}/include)
    target_compile_definitions(wirehair_objects PRIVATE WIREHAIR_BUILDING=1)
    target_link_libraries(wirehair_objects PRIVATE Threads::Threads)
    wirehair_apply_project_options(wirehair_objects)

    add_library(wirehair ${WIREHAIR_PRIMARY_VARIANT}
        $<TARGET_OBJECTS:wirehair_objects>)
    wirehair_configure_library(wirehair ${WIREHAIR_PRIMARY_VARIANT})
    if(BUILD_SHARED_LIBS)
        add_library(wirehair_static STATIC $<TARGET_OBJECTS:wirehair_objects>)
        set_target_properties(wirehair_static PROPERTIES OUTPUT_NAME wirehair)
        wirehair_configure_library(wirehair_static STATIC)
        list(APPEND WIREHAIR_INSTALL_TARGETS wirehair_static)
    else()
        add_library(wirehair_shared SHARED $<TARGET_OBJECTS:wirehair_objects>)
        set_target_properties(wirehair_shared PROPERTIES OUTPUT_NAME wirehair)
        wirehair_configure_library(wirehair_shared SHARED)
        list(APPEND WIREHAIR_INSTALL_TARGETS wirehair_shared)
    endif()
else()
    add_library(wirehair ${WIREHAIR_PRIMARY_VARIANT} ${LIB_SOURCE_FILES})
    wirehair_configure_library(wirehair ${WIREHAIR_PRIMARY_VARIANT})
    if(WIREHAIR_BUILD_BOTH)
        # Windows static and DLL objects require different export definitions.
        if(BUILD_SHARED_LIBS)
            add_library(wirehair_static STATIC ${LIB_SOURCE_FILES})
            set_target_properties(wirehair_static PROPERTIES
                OUTPUT_NAME wirehair_static)
            wirehair_configure_library(wirehair_static STATIC)
            list(APPEND WIREHAIR_INSTALL_TARGETS wirehair_static)
        else()
            add_library(wirehair_shared SHARED ${LIB_SOURCE_FILES})
            set_target_properties(wirehair_shared PROPERTIES
                OUTPUT_NAME wirehair_shared
                # Keep the import library distinct from the primary static
                # archive, while giving the DLL the standard name searched by
                # the installed ctypes binding.
                RUNTIME_OUTPUT_NAME wirehair)
            wirehair_configure_library(wirehair_shared SHARED)
            list(APPEND WIREHAIR_INSTALL_TARGETS wirehair_shared)
        endif()
    endif()
endif()
add_library(wirehair::wirehair ALIAS wirehair)

if (BUILD_TESTS)
    enable_testing()
    find_package(OpenMP QUIET)

    add_executable(unit_test ${UNIT_TEST_SOURCE_FILES})
    wirehair_apply_project_options(unit_test)
    target_link_libraries(unit_test wirehair)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(unit_test OpenMP::OpenMP_CXX)
    endif()
    add_test(NAME unit_test COMMAND unit_test)
    set_tests_properties(unit_test PROPERTIES
        ENVIRONMENT "WIREHAIR_UNIT_SKIP_BENCH=1;WIREHAIR_UNIT_MIN_N=2;WIREHAIR_UNIT_MAX_N=256;WIREHAIR_UNIT_BLOCK_BYTES_MIN=1;WIREHAIR_UNIT_BLOCK_BYTES_MAX=32")

    add_test(NAME unit_test_terminal_failure
        COMMAND ${CMAKE_COMMAND}
            -DBINARY=$<TARGET_FILE:unit_test>
            -DEXPECT_EXIT=4
            "-DEXPECT_PATTERN=Terminal decoder exhaustion"
            -P ${PROJECT_SOURCE_DIR}/test/ExpectFailure.cmake)
    set_tests_properties(unit_test_terminal_failure PROPERTIES
        ENVIRONMENT "WIREHAIR_UNIT_SKIP_BENCH=1;WIREHAIR_UNIT_MIN_N=2;WIREHAIR_UNIT_MAX_N=2;WIREHAIR_UNIT_BLOCK_BYTES_MIN=1;WIREHAIR_UNIT_BLOCK_BYTES_MAX=1;WIREHAIR_UNIT_INJECT_TERMINAL=1"
        TIMEOUT 15)

    # Fault injection is confined to this non-installed test library.
    add_library(wirehair_test_support STATIC ${LIB_SOURCE_FILES})
    target_include_directories(wirehair_test_support PUBLIC ${PROJECT_SOURCE_DIR}/include)
    target_compile_definitions(wirehair_test_support
        PRIVATE WIREHAIR_TESTING=1 WIREHAIR_BUILDING=1
        PUBLIC WIREHAIR_STATIC=1)
    target_link_libraries(wirehair_test_support PRIVATE Threads::Threads)
    set_target_properties(wirehair_test_support PROPERTIES
        POSITION_INDEPENDENT_CODE ${WIREHAIR_STATIC_PIC})
    wirehair_apply_project_options(wirehair_test_support)

    add_executable(v2_borrowed_facade_fault_test
        test/V2BorrowedFacadeFaultTest.cpp)
    wirehair_apply_project_options(v2_borrowed_facade_fault_test)
    target_link_libraries(v2_borrowed_facade_fault_test
        wirehair_test_support Threads::Threads)
    add_test(NAME v2_borrowed_facade_fault_test
        COMMAND v2_borrowed_facade_fault_test)
    set_tests_properties(v2_borrowed_facade_fault_test PROPERTIES
        LABELS "v2;borrowed-source;allocation;transactional"
        TIMEOUT 120)

    add_executable(legacy_core_test test/LegacyCoreTest.cpp)
    wirehair_apply_project_options(legacy_core_test)
    target_link_libraries(legacy_core_test wirehair_test_support Threads::Threads)
    add_test(NAME legacy_core_test COMMAND legacy_core_test)
    add_test(NAME legacy_core_preinit_test COMMAND legacy_core_test --preinit)
    add_test(NAME legacy_core_init_failure_test COMMAND legacy_core_test --init-failure)
    add_test(NAME legacy_core_cold_start_test COMMAND legacy_core_test --cold-start)
    add_test(NAME legacy_duplicate_packet_test COMMAND legacy_core_test --duplicates)
    set_tests_properties(legacy_duplicate_packet_test PROPERTIES
        LABELS "decoder;duplicates;networking;e2e"
        TIMEOUT 120)
    add_test(NAME legacy_encoder_detach_test COMMAND legacy_core_test --detach)
    set_tests_properties(legacy_encoder_detach_test PROPERTIES
        LABELS "encoder;lifetime;memory;e2e"
        TIMEOUT 120)

    add_executable(peel_avalanche_test test/PeelAvalancheTest.cpp)
    wirehair_apply_project_options(peel_avalanche_test)
    target_link_libraries(peel_avalanche_test
        wirehair_test_support Threads::Threads)
    add_test(NAME peel_avalanche_small_stack
        COMMAND peel_avalanche_test --small-stack)
    set_tests_properties(peel_avalanche_small_stack PROPERTIES
        LABELS "decoder;peeling;memory-safety;e2e"
        TIMEOUT 120
        PROCESSORS 1)

    add_executable(gf256_inplace_test test/GF256InPlaceTest.cpp)
    wirehair_apply_project_options(gf256_inplace_test)
    target_link_libraries(gf256_inplace_test wirehair_test_support)
    add_test(NAME gf256_inplace_test COMMAND gf256_inplace_test)
    set_tests_properties(gf256_inplace_test PROPERTIES
        LABELS "gf256;simd;undefined-behavior")

    add_executable(heavy_window_test test/HeavyWindowTest.cpp)
    wirehair_apply_project_options(heavy_window_test)
    add_test(NAME heavy_window_test COMMAND heavy_window_test)
    set_tests_properties(heavy_window_test PROPERTIES
        LABELS "heavy;portability;property")

    add_executable(portability_roundtrip_test
        test/PortabilityRoundTripTest.cpp)
    wirehair_apply_project_options(portability_roundtrip_test)
    target_link_libraries(portability_roundtrip_test wirehair)
    add_test(NAME portability_roundtrip_test
        COMMAND portability_roundtrip_test)
    set_tests_properties(portability_roundtrip_test PROPERTIES
        LABELS "portability;golden;e2e"
        TIMEOUT 60)

    add_executable(allocation_overflow_test test/AllocationOverflowTest.cpp)
    wirehair_apply_project_options(allocation_overflow_test)
    target_compile_definitions(allocation_overflow_test
        PRIVATE WIREHAIR_TESTING=1)
    target_link_libraries(allocation_overflow_test
        wirehair_test_support Threads::Threads)
    add_test(NAME allocation_overflow_test COMMAND allocation_overflow_test)
    set_tests_properties(allocation_overflow_test PROPERTIES
        LABELS "allocation;memory-safety;32-bit"
        ENVIRONMENT "ASAN_OPTIONS=allocator_may_return_null=1:detect_leaks=1:halt_on_error=1:strict_string_checks=1;UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1"
        TIMEOUT 60)

    add_executable(legacy_c_consumer test/LegacyCConsumer.c)
    set_target_properties(legacy_c_consumer PROPERTIES LINKER_LANGUAGE CXX)
    wirehair_apply_project_options(legacy_c_consumer)
    target_link_libraries(legacy_c_consumer wirehair)
    add_test(NAME legacy_c_consumer_test COMMAND legacy_c_consumer)
    add_executable(legacy_cpp_consumer test/LegacyCppConsumer.cpp)
    wirehair_apply_project_options(legacy_cpp_consumer)
    target_link_libraries(legacy_cpp_consumer wirehair)
    add_test(NAME legacy_cpp_consumer_test COMMAND legacy_cpp_consumer)
    if(TARGET wirehair_shared)
        add_executable(legacy_c_consumer_shared test/LegacyCConsumer.c)
        set_target_properties(legacy_c_consumer_shared PROPERTIES LINKER_LANGUAGE CXX)
        wirehair_apply_project_options(legacy_c_consumer_shared)
        target_link_libraries(legacy_c_consumer_shared wirehair_shared)
        add_test(NAME legacy_c_consumer_shared_test COMMAND legacy_c_consumer_shared)
    endif()

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(WIREHAIR_ELF_SHARED_TARGET "")
        if(WIREHAIR_PRIMARY_VARIANT STREQUAL "SHARED")
            set(WIREHAIR_ELF_SHARED_TARGET wirehair)
        elseif(TARGET wirehair_shared)
            set(WIREHAIR_ELF_SHARED_TARGET wirehair_shared)
        endif()
        if(WIREHAIR_ELF_SHARED_TARGET)
            add_test(
                NAME elf_public_exports
                COMMAND ${CMAKE_COMMAND}
                    -DLIBRARY=$<TARGET_FILE:${WIREHAIR_ELF_SHARED_TARGET}>
                    -DMANIFEST=${WIREHAIR_ELF_EXPORT_MAP}
                    -DNM=${CMAKE_NM}
                    -P ${PROJECT_SOURCE_DIR}/test/cmake/CheckElfExports.cmake)
            set_tests_properties(elf_public_exports PROPERTIES
                LABELS "abi;elf;exports"
                TIMEOUT 30)
            add_test(
                NAME elf_public_exports_negative
                COMMAND ${CMAKE_COMMAND}
                    -DCHECKER=${PROJECT_SOURCE_DIR}/test/cmake/CheckElfExports.cmake
                    -DLIBRARY=$<TARGET_FILE:${WIREHAIR_ELF_SHARED_TARGET}>
                    -DMANIFEST=${WIREHAIR_ELF_EXPORT_MAP}
                    -DNM=${CMAKE_NM}
                    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/elf-export-negative
                    -P ${PROJECT_SOURCE_DIR}/test/cmake/CheckElfExportsSelfTest.cmake)
            set_tests_properties(elf_public_exports_negative PROPERTIES
                DEPENDS elf_public_exports
                LABELS "abi;elf;exports;negative"
                TIMEOUT 30)
        endif()
    endif()

    add_executable(large_message_test ${LARGE_MESSAGE_TEST_SOURCE_FILES})
    wirehair_apply_project_options(large_message_test)
    target_link_libraries(large_message_test wirehair)

    set(WIREHAIR_LARGE_FAST_ENV
        "WIREHAIR_LARGE_N=256;WIREHAIR_LARGE_BLOCK_BYTES=257;WIREHAIR_LARGE_FINAL_BYTES=113;WIREHAIR_LARGE_LOSS_COUNT=5;WIREHAIR_LARGE_MAX_REPAIR=64;WIREHAIR_LARGE_SEED=0x8af09d31e7c4526b;WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=1048576;WIREHAIR_LARGE_MAX_MILLISECONDS=30000")
    add_test(NAME large_message_fast COMMAND large_message_test)
    set_tests_properties(large_message_fast PROPERTIES
        ENVIRONMENT "${WIREHAIR_LARGE_FAST_ENV}"
        TIMEOUT 45
        LABELS "e2e;large-message;fast"
        PROCESSORS 1)

    add_test(NAME large_message_terminal_failure
        COMMAND ${CMAKE_COMMAND}
            -DBINARY=$<TARGET_FILE:large_message_test>
            -DEXPECT_EXIT=1
            "-DEXPECT_PATTERN=Injected terminal decode result"
            -P ${PROJECT_SOURCE_DIR}/test/ExpectFailure.cmake)
    set_tests_properties(large_message_terminal_failure PROPERTIES
        ENVIRONMENT "${WIREHAIR_LARGE_FAST_ENV};WIREHAIR_LARGE_INJECT_TERMINAL=1"
        TIMEOUT 15
        LABELS "e2e;large-message;failure"
        PROCESSORS 1)

    set(WIREHAIR_LARGE_REPAIR_CAP_ENV
        "WIREHAIR_LARGE_N=256;WIREHAIR_LARGE_BLOCK_BYTES=257;WIREHAIR_LARGE_FINAL_BYTES=113;WIREHAIR_LARGE_LOSS_COUNT=8;WIREHAIR_LARGE_MAX_REPAIR=1;WIREHAIR_LARGE_SEED=0x8af09d31e7c4526b;WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=1048576;WIREHAIR_LARGE_MAX_MILLISECONDS=30000")
    add_test(NAME large_message_repair_cap_failure
        COMMAND ${CMAKE_COMMAND}
            -DBINARY=$<TARGET_FILE:large_message_test>
            -DEXPECT_EXIT=1
            "-DEXPECT_PATTERN=repair blocks \\(cap reached\\)"
            -P ${PROJECT_SOURCE_DIR}/test/ExpectFailure.cmake)
    set_tests_properties(large_message_repair_cap_failure PROPERTIES
        ENVIRONMENT "${WIREHAIR_LARGE_REPAIR_CAP_ENV}"
        TIMEOUT 15
        LABELS "e2e;large-message;failure"
        PROCESSORS 1)

    add_test(NAME large_message_timeout_failure
        COMMAND ${CMAKE_COMMAND}
            -DBINARY=$<TARGET_FILE:large_message_test>
            -DEXPECT_EXIT=3
            "-DEXPECT_PATTERN=Injected timeout failure"
            -P ${PROJECT_SOURCE_DIR}/test/ExpectFailure.cmake)
    set_tests_properties(large_message_timeout_failure PROPERTIES
        ENVIRONMENT "${WIREHAIR_LARGE_FAST_ENV};WIREHAIR_LARGE_INJECT_TIMEOUT=1"
        TIMEOUT 15
        LABELS "e2e;large-message;failure"
        PROCESSORS 1)

    if(UNIX)
        add_test(NAME large_message_runner_contract
            COMMAND bash ${PROJECT_SOURCE_DIR}/test/test_run_large_message_profiles.sh)
        set_tests_properties(large_message_runner_contract PROPERTIES
            TIMEOUT 30
            LABELS "large-message;runner;contract")
    endif()

    if(WIREHAIR_ENABLE_SCHEDULED_TESTS)
        add_test(NAME large_message_n64000_byte COMMAND large_message_test)
        set_tests_properties(large_message_n64000_byte PROPERTIES
            ENVIRONMENT "WIREHAIR_LARGE_N=64000;WIREHAIR_LARGE_BLOCK_BYTES=1;WIREHAIR_LARGE_FINAL_BYTES=1;WIREHAIR_LARGE_LOSS_COUNT=17;WIREHAIR_LARGE_MAX_REPAIR=128;WIREHAIR_LARGE_SEED=0x640001;WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=65536;WIREHAIR_LARGE_MAX_MILLISECONDS=300000"
            TIMEOUT 330
            RUN_SERIAL TRUE
            PROCESSORS 1
            LABELS "e2e;large-message;scheduled;high-n;one-byte")

        add_test(NAME large_message_packet_partial COMMAND large_message_test)
        set_tests_properties(large_message_packet_partial PROPERTIES
            ENVIRONMENT "WIREHAIR_LARGE_N=4096;WIREHAIR_LARGE_BLOCK_BYTES=1280;WIREHAIR_LARGE_FINAL_BYTES=733;WIREHAIR_LARGE_LOSS_COUNT=11;WIREHAIR_LARGE_MAX_REPAIR=96;WIREHAIR_LARGE_SEED=0x40961280;WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=8388608;WIREHAIR_LARGE_MAX_MILLISECONDS=180000"
            TIMEOUT 210
            RUN_SERIAL TRUE
            PROCESSORS 1
            LABELS "e2e;large-message;scheduled;packet;partial-final")

        add_test(NAME large_message_large_block_partial COMMAND large_message_test)
        set_tests_properties(large_message_large_block_partial PROPERTIES
            ENVIRONMENT "WIREHAIR_LARGE_N=32;WIREHAIR_LARGE_BLOCK_BYTES=1048576;WIREHAIR_LARGE_FINAL_BYTES=1048559;WIREHAIR_LARGE_LOSS_COUNT=3;WIREHAIR_LARGE_MAX_REPAIR=64;WIREHAIR_LARGE_SEED=0x321048576;WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=33554432;WIREHAIR_LARGE_MAX_MILLISECONDS=180000"
            TIMEOUT 210
            RUN_SERIAL TRUE
            PROCESSORS 1
            LABELS "e2e;large-message;scheduled;large-block;partial-final;memory")
    endif()
endif()

if(WIREHAIR_BUILD_TOOLS)
    set(WIREHAIR_TOOL_EXCLUDE "")
else()
    set(WIREHAIR_TOOL_EXCLUDE EXCLUDE_FROM_ALL)
endif()

# These targets remain directly buildable even when excluded from `all`.
find_package(OpenMP QUIET)
set(WIREHAIR_INTERNAL_LIBRARY wirehair)
if(WIREHAIR_PRIMARY_VARIANT STREQUAL "SHARED")
    if(TARGET wirehair_static)
        set(WIREHAIR_INTERNAL_LIBRARY wirehair_static)
    else()
        # Shared libraries expose the public C ABI only.  V2 and offline
        # generators use internal C++ codec types, so give those targets an
        # explicit static implementation rather than relying on leaked DSO
        # symbols.  It stays out of routine builds unless an internal target
        # needs it.
        add_library(wirehair_internal_support STATIC EXCLUDE_FROM_ALL
            ${LIB_SOURCE_FILES})
        wirehair_configure_library(wirehair_internal_support STATIC)
        set(WIREHAIR_INTERNAL_LIBRARY wirehair_internal_support)
    endif()
endif()

# Coverage-guided public-facade fuzzing needs feedback from the production
# implementation itself, not only from the harness translation unit.  Build a
# non-installed, production-equivalent static copy with sanitizer coverage;
# do not enable the fault-injection hooks used by wirehair_test_support.
set(WIREHAIR_PUBLIC_FUZZ_LIBRARY wirehair)
if(BUILD_TESTS AND BUILD_CODEC_V2 AND WIREHAIR_ENABLE_LIBFUZZER)
    add_library(wirehair_public_fuzz_support STATIC EXCLUDE_FROM_ALL
        ${LIB_SOURCE_FILES})
    wirehair_configure_library(wirehair_public_fuzz_support STATIC)
    target_compile_options(wirehair_public_fuzz_support PRIVATE
        -fsanitize=fuzzer-no-link)
    set(WIREHAIR_PUBLIC_FUZZ_LIBRARY wirehair_public_fuzz_support)
endif()
add_executable(gen_small_dseeds ${WIREHAIR_TOOL_EXCLUDE} ${GEN_SMALL_DSEEDS})
add_executable(gen_peel_seeds ${WIREHAIR_TOOL_EXCLUDE} ${GEN_PEEL_SEEDS})
add_executable(gen_most_dseeds ${WIREHAIR_TOOL_EXCLUDE} ${GEN_MOST_DSEEDS})
add_executable(gen_dcounts ${WIREHAIR_TOOL_EXCLUDE} ${GEN_DCOUNTS})
foreach(_wirehair_generator IN ITEMS
        gen_small_dseeds gen_peel_seeds gen_most_dseeds gen_dcounts)
    target_link_libraries(${_wirehair_generator} PRIVATE
        ${WIREHAIR_INTERNAL_LIBRARY})
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${_wirehair_generator} PRIVATE OpenMP::OpenMP_CXX)
    endif()
    wirehair_apply_project_options(${_wirehair_generator})
endforeach()
add_executable(gen_tables ${WIREHAIR_TOOL_EXCLUDE} ${GEN_TABLES})
target_link_libraries(gen_tables PRIVATE Threads::Threads)
wirehair_apply_project_options(gen_tables)

if (BUILD_CODEC_V2)
    add_subdirectory(codec)
endif()

# The header is installed by the include/ DIRECTORY rule below; the target
# sets no PUBLIC_HEADER property, so a PUBLIC_HEADER clause here would be a
# silent no-op.
install(TARGETS ${WIREHAIR_INSTALL_TARGETS}
    EXPORT wirehairTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# The ctypes binding cannot use a static archive.  Install it only when this
# package also contains a shared library that it can load at runtime.
if(BUILD_SHARED_LIBS OR WIREHAIR_BUILD_BOTH)
    install(FILES python/whirehair.py DESTINATION python)
    install(DIRECTORY python/wirehair DESTINATION python
        FILES_MATCHING
        PATTERN "*.py"
        PATTERN "__pycache__" EXCLUDE)
endif()

set(WIREHAIR_PACKAGE_INSTALL_DIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/wirehair")

# Generate a relocatable pkg-config file beside the installed library.  The
# path from <libdir>/pkgconfig back to the prefix depends on how deeply a
# caller nests CMAKE_INSTALL_LIBDIR, so calculate it instead of assuming
# "lib/pkgconfig".
set(WIREHAIR_PKGCONFIG_INSTALL_DIR
    "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
if(IS_ABSOLUTE "${WIREHAIR_PKGCONFIG_INSTALL_DIR}")
    set(WIREHAIR_PC_PREFIX "${CMAKE_INSTALL_PREFIX}")
    set(WIREHAIR_PC_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
else()
    # file(RELATIVE_PATH) requires native absolute paths.  A Unix-rooted
    # synthetic path is not absolute on Windows (it has no drive), so anchor
    # both operands beneath this configure tree on every host.  The anchor
    # need not exist; only the lexical relative path is used in wirehair.pc.
    set(_wirehair_pc_prefix_anchor
        "${CMAKE_CURRENT_BINARY_DIR}/wirehair-pc-prefix-anchor")
    file(RELATIVE_PATH WIREHAIR_PC_PREFIX_RELATIVE
        "${_wirehair_pc_prefix_anchor}/${WIREHAIR_PKGCONFIG_INSTALL_DIR}"
        "${_wirehair_pc_prefix_anchor}")
    unset(_wirehair_pc_prefix_anchor)
    set(WIREHAIR_PC_PREFIX
        "\${pcfiledir}/${WIREHAIR_PC_PREFIX_RELATIVE}")
    set(WIREHAIR_PC_LIBDIR "\${pcfiledir}/..")
endif()
if(IS_ABSOLUTE "${CMAKE_INSTALL_INCLUDEDIR}")
    set(WIREHAIR_PC_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}")
else()
    set(WIREHAIR_PC_INCLUDEDIR
        "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
endif()

# Static consumers compiling as C must pull in the C++ runtime, threading, and
# math libraries used by the implementation.  Keep those dependencies private
# for shared/dual installs; a static-only install exposes them in Libs because
# even plain `--libs` can only select the archive in that package.
list(FIND CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++"
    WIREHAIR_PC_LIBSTDCXX_INDEX)
list(FIND CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "c++"
    WIREHAIR_PC_LIBCXX_INDEX)
if(NOT WIREHAIR_PC_LIBSTDCXX_INDEX EQUAL -1)
    set(WIREHAIR_PC_RUNTIME_LIBS "-lstdc++ -pthread -lm")
elseif(NOT WIREHAIR_PC_LIBCXX_INDEX EQUAL -1)
    set(WIREHAIR_PC_RUNTIME_LIBS "-lc++ -pthread -lm")
elseif(APPLE)
    set(WIREHAIR_PC_RUNTIME_LIBS "-lc++ -pthread -lm")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|Intel")
    set(WIREHAIR_PC_RUNTIME_LIBS "-lstdc++ -pthread -lm")
elseif(UNIX OR MINGW)
    message(FATAL_ERROR
        "Cannot determine the C++ runtime needed by static pkg-config users; "
        "implicit libraries are '${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}'")
else()
    set(WIREHAIR_PC_RUNTIME_LIBS "")
endif()
if(NOT BUILD_SHARED_LIBS AND NOT WIREHAIR_BUILD_BOTH)
    set(WIREHAIR_PC_PUBLIC_RUNTIME_LIBS "${WIREHAIR_PC_RUNTIME_LIBS}")
    set(WIREHAIR_PC_PRIVATE_RUNTIME_LIBS "")
else()
    set(WIREHAIR_PC_PUBLIC_RUNTIME_LIBS "")
    set(WIREHAIR_PC_PRIVATE_RUNTIME_LIBS "${WIREHAIR_PC_RUNTIME_LIBS}")
endif()
configure_file(
    cmake/wirehair.pc.in
    "${CMAKE_CURRENT_BINARY_DIR}/wirehair.pc"
    @ONLY)
configure_package_config_file(
    cmake/wirehairConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/wirehairConfig.cmake"
    INSTALL_DESTINATION "${WIREHAIR_PACKAGE_INSTALL_DIR}")
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/wirehairConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
install(EXPORT wirehairTargets
    FILE wirehairTargets.cmake
    NAMESPACE wirehair::
    DESTINATION "${WIREHAIR_PACKAGE_INSTALL_DIR}")
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/wirehairConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/wirehairConfigVersion.cmake"
    DESTINATION "${WIREHAIR_PACKAGE_INSTALL_DIR}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/wirehair.pc"
    DESTINATION "${WIREHAIR_PKGCONFIG_INSTALL_DIR}")
install(FILES LICENSE
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/wirehair")
install(FILES
    LEGACY_WIRE_PROFILES.md
    V2_WIRE_PROFILE.md
    DESTINATION "${CMAKE_INSTALL_DATADIR}/doc/wirehair")

if(BUILD_TESTS)
    find_package(Python3 3.8 COMPONENTS Interpreter REQUIRED)
    if(WIREHAIR_ENABLE_TRI_BMQ4_SCHUR_FALSIFIER)
        add_executable(wirehair_wh2_tri_bmq4_schur_falsifier
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2TriBmq4SchurFalsifier.cpp)
        target_link_libraries(
            wirehair_wh2_tri_bmq4_schur_falsifier PRIVATE
            wirehair_v2_policy)
        target_compile_definitions(
            wirehair_wh2_tri_bmq4_schur_falsifier PRIVATE
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\"")
        wirehair_apply_project_options(
            wirehair_wh2_tri_bmq4_schur_falsifier)
        add_test(
            NAME wirehair_wh2_tri_bmq4_schur_falsifier_selftest
            COMMAND wirehair_wh2_tri_bmq4_schur_falsifier --selftest)
        set_tests_properties(
            wirehair_wh2_tri_bmq4_schur_falsifier_selftest PROPERTIES
                TIMEOUT 60
                LABELS "wh2;experimental;structure;falsifier")
    endif()
    add_executable(wirehair_wh2_frozen_trace_test
        bench/Wh2FrozenTrace.cpp
        bench/Wh2FrozenTrace.h
        bench/Wh2FrozenTraceTest.cpp)
    wirehair_apply_project_options(wirehair_wh2_frozen_trace_test)
    add_test(
        NAME wirehair_wh2_frozen_trace_test
        COMMAND wirehair_wh2_frozen_trace_test)
    set_tests_properties(wirehair_wh2_frozen_trace_test PROPERTIES
        TIMEOUT 120
        LABELS "wh2;benchmark;contract;native;trace")
    add_executable(wirehair_wh2_frozen_timing_test
        bench/Wh2FrozenTrace.cpp
        bench/Wh2FrozenTrace.h
        bench/Wh2FrozenTiming.cpp
        bench/Wh2FrozenTiming.h
        bench/Wh2FrozenTimingTest.cpp)
    wirehair_apply_project_options(wirehair_wh2_frozen_timing_test)
    add_test(
        NAME wirehair_wh2_frozen_timing_test
        COMMAND wirehair_wh2_frozen_timing_test)
    set_tests_properties(wirehair_wh2_frozen_timing_test PROPERTIES
        TIMEOUT 120
        LABELS "wh2;benchmark;contract;native;timing")
    add_executable(wirehair_wh2_native_panel_test
        bench/Wh2NativePanel.cpp
        bench/Wh2NativePanel.h
        bench/Wh2NativePanelTest.cpp)
    wirehair_apply_project_options(wirehair_wh2_native_panel_test)
    add_test(
        NAME wirehair_wh2_native_panel_test
        COMMAND wirehair_wh2_native_panel_test)
    set_tests_properties(wirehair_wh2_native_panel_test PROPERTIES
        TIMEOUT 120
        LABELS "wh2;benchmark;contract;native;timing;affinity")
    add_executable(wirehair_wh2_public_borrowed_current_vs_wh1_r1
        bench/Wh2PublicBorrowedCurrentVsWh1R1.cpp
        bench/Wh2FrozenTrace.cpp
        bench/Wh2FrozenTrace.h
        bench/Wh2NativePanel.cpp
        bench/Wh2NativePanel.h
        bench/Wh2PublicBorrowedTargetIdentity.cpp
        bench/Wh2PublicBorrowedTargetIdentity.h
        bench/Wh2RdpruTargetIdentityV2.cpp
        bench/Wh2RdpruTargetIdentityV2.h)
    get_filename_component(
        WIREHAIR_WH2_CXX_COMPILER_REALPATH
        "${CMAKE_CXX_COMPILER}" REALPATH)
    file(SHA256 "${WIREHAIR_WH2_CXX_COMPILER_REALPATH}"
        WIREHAIR_WH2_CXX_COMPILER_SHA256)
    target_link_libraries(
        wirehair_wh2_public_borrowed_current_vs_wh1_r1 PRIVATE
        wirehair
        ${CMAKE_DL_LIBS})
    target_compile_definitions(
        wirehair_wh2_public_borrowed_current_vs_wh1_r1 PRIVATE
        "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\""
        "WIREHAIR_WH2_CXX_COMPILER_PATH=\"${WIREHAIR_WH2_CXX_COMPILER_REALPATH}\""
        "WIREHAIR_WH2_CXX_COMPILER_SHA256=\"${WIREHAIR_WH2_CXX_COMPILER_SHA256}\""
        "WIREHAIR_WH2_CMAKE_CXX_COMPILER_ID=\"${CMAKE_CXX_COMPILER_ID}\""
        "WIREHAIR_WH2_CMAKE_CXX_COMPILER_VERSION=\"${CMAKE_CXX_COMPILER_VERSION}\"")
    wirehair_apply_project_options(
        wirehair_wh2_public_borrowed_current_vs_wh1_r1)
    add_test(
        NAME wirehair_wh2_public_borrowed_current_vs_wh1_r1_selftest
        COMMAND wirehair_wh2_public_borrowed_current_vs_wh1_r1 --selftest)
    set_tests_properties(
        wirehair_wh2_public_borrowed_current_vs_wh1_r1_selftest PROPERTIES
            TIMEOUT 120
            LABELS
                "wh2;benchmark;contract;public-api;encoder;facade;timing;r1")
    add_test(
        NAME wirehair_wh2_public_borrowed_current_vs_wh1_r1_python_selftest
        COMMAND ${Python3_EXECUTABLE} -I -B -S
            ${PROJECT_SOURCE_DIR}/bench/Wh2PublicBorrowedCurrentVsWh1R1.py
            --selftest)
    set_tests_properties(
        wirehair_wh2_public_borrowed_current_vs_wh1_r1_python_selftest
        PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 60
            LABELS
                "wh2;benchmark;contract;public-api;encoder;facade;timing;r1")
    add_test(
        NAME wirehair_wh2_public_borrowed_current_vs_wh1_r1_python_test
        COMMAND ${Python3_EXECUTABLE} -I -B -S
            ${PROJECT_SOURCE_DIR}/bench/test_Wh2PublicBorrowedCurrentVsWh1R1.py)
    set_tests_properties(
        wirehair_wh2_public_borrowed_current_vs_wh1_r1_python_test
        PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS
                "wh2;benchmark;contract;public-api;encoder;facade;timing;r1")
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
       NOT CMAKE_CROSSCOMPILING AND BUILD_SHARED_LIBS)
        if(NOT CMAKE_NM)
            message(FATAL_ERROR
                "native Linux shared public-borrowed r1 contract requires CMAKE_NM")
        endif()
        add_test(
            NAME wirehair_wh2_public_borrowed_current_vs_wh1_r1_symbol_audit
            COMMAND ${CMAKE_COMMAND}
                "-DNM=${CMAKE_NM}"
                "-DPUBLIC_BORROWED_R1=$<TARGET_FILE:wirehair_wh2_public_borrowed_current_vs_wh1_r1>"
                -P ${PROJECT_SOURCE_DIR}/cmake/Wh2PublicBorrowedCurrentVsWh1R1SymbolAudit.cmake)
        set_tests_properties(
            wirehair_wh2_public_borrowed_current_vs_wh1_r1_symbol_audit
            PROPERTIES
                TIMEOUT 30
                LABELS
                    "wh2;benchmark;contract;public-api;encoder;facade;timing;r1;symbols")
        # This oracle pins the preregistered CPU and captures its exact
        # capabilities. Generic CI hosts need not have CPU 120 or GCC 13.
        # The real controller always runs the oracle during preflight,
        # independently of whether it is also registered with CTest.
        if(WIREHAIR_ENABLE_PUBLIC_BORROWED_R1_NATIVE_CONFIG)
            add_test(
                NAME wirehair_wh2_public_borrowed_current_vs_wh1_r1_native_config
                COMMAND ${Python3_EXECUTABLE} -I -B -S
                    ${PROJECT_SOURCE_DIR}/bench/Wh2PublicBorrowedCurrentVsWh1R1.py
                    --native-config-selftest
                    --worker $<TARGET_FILE:wirehair_wh2_public_borrowed_current_vs_wh1_r1>
                    --library $<TARGET_FILE:wirehair>
                    --compiler ${WIREHAIR_WH2_CXX_COMPILER_REALPATH}
                    --expected-source-commit ${WIREHAIR_WH2_SOURCE_GIT_COMMIT})
            set_tests_properties(
                wirehair_wh2_public_borrowed_current_vs_wh1_r1_native_config
                PROPERTIES
                    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                    ENVIRONMENT "PYTHONWARNINGS=error"
                    RUN_SERIAL TRUE
                    TIMEOUT 360
                    LABELS
                        "wh2;benchmark;contract;public-api;encoder;facade;r1;native-config")
        endif()
    endif()
    # Independent A/A calibration. The main TU includes the immutable r1
    # helper implementation under a renamed entrypoint; it never runs r1.
    add_executable(wirehair_wh2_public_aa_calibration_r0
        bench/Wh2PublicAaCalibrationR0.cpp
        bench/Wh2FrozenTrace.cpp
        bench/Wh2FrozenTrace.h
        bench/Wh2NativePanel.cpp
        bench/Wh2NativePanel.h
        bench/Wh2PublicBorrowedTargetIdentity.cpp
        bench/Wh2PublicBorrowedTargetIdentity.h
        bench/Wh2RdpruTargetIdentityV2.cpp
        bench/Wh2RdpruTargetIdentityV2.h)
    target_link_libraries(wirehair_wh2_public_aa_calibration_r0 PRIVATE
        wirehair ${CMAKE_DL_LIBS})
    get_target_property(_wh2_aa_calibration_definitions
        wirehair_wh2_public_borrowed_current_vs_wh1_r1 COMPILE_DEFINITIONS)
    target_compile_definitions(wirehair_wh2_public_aa_calibration_r0 PRIVATE
        ${_wh2_aa_calibration_definitions})
    wirehair_apply_project_options(wirehair_wh2_public_aa_calibration_r0)
    add_test(NAME wirehair_wh2_public_aa_calibration_r0_selftest
        COMMAND wirehair_wh2_public_aa_calibration_r0 --selftest)
    add_test(NAME wirehair_wh2_public_aa_calibration_r0_python_selftest
        COMMAND ${Python3_EXECUTABLE} -I -B -S
            ${PROJECT_SOURCE_DIR}/bench/Wh2PublicAaCalibrationR0.py --selftest)
    add_test(NAME wirehair_wh2_public_aa_calibration_r0_python_test
        COMMAND ${Python3_EXECUTABLE} -I -B -S
            ${PROJECT_SOURCE_DIR}/bench/test_Wh2PublicAaCalibrationR0.py)
    set_tests_properties(
        wirehair_wh2_public_aa_calibration_r0_selftest
        wirehair_wh2_public_aa_calibration_r0_python_selftest
        wirehair_wh2_public_aa_calibration_r0_python_test
        PROPERTIES TIMEOUT 120 WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        ENVIRONMENT "PYTHONWARNINGS=error"
        LABELS "wh2;benchmark;contract;public-api;aa-calibration")
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
       NOT CMAKE_CROSSCOMPILING AND BUILD_SHARED_LIBS)
        add_test(NAME wirehair_wh2_public_aa_calibration_r0_symbol_audit
            COMMAND ${CMAKE_COMMAND}
                "-DNM=${CMAKE_NM}"
                "-DPUBLIC_BORROWED_R1=$<TARGET_FILE:wirehair_wh2_public_aa_calibration_r0>"
                -P ${PROJECT_SOURCE_DIR}/cmake/Wh2PublicBorrowedCurrentVsWh1R1SymbolAudit.cmake)
        set_tests_properties(wirehair_wh2_public_aa_calibration_r0_symbol_audit
            PROPERTIES TIMEOUT 30
            LABELS "wh2;benchmark;contract;public-api;aa-calibration;symbols")
    endif()
    # Small page-layout A/A mechanism screen; public APIs only.
    add_executable(wirehair_wh2_public_page_aa_r0
        bench/Wh2PublicPageAaR0.cpp
        bench/Wh2FrozenTrace.cpp
        bench/Wh2NativePanel.cpp
        bench/Wh2PublicBorrowedTargetIdentity.cpp
        bench/Wh2RdpruTargetIdentityV2.cpp)
    target_link_libraries(wirehair_wh2_public_page_aa_r0 PRIVATE
        wirehair ${CMAKE_DL_LIBS})
    target_compile_definitions(wirehair_wh2_public_page_aa_r0 PRIVATE
        ${_wh2_aa_calibration_definitions})
    wirehair_apply_project_options(wirehair_wh2_public_page_aa_r0)
    add_test(NAME wirehair_wh2_public_page_aa_r0_selftest
        COMMAND wirehair_wh2_public_page_aa_r0 --selftest)
    add_test(NAME wirehair_wh2_public_page_aa_r0_python_selftest
        COMMAND ${Python3_EXECUTABLE} -I -B -S
            ${PROJECT_SOURCE_DIR}/bench/Wh2PublicPageAaR0.py --selftest)
    add_test(NAME wirehair_wh2_public_page_aa_r0_python_test
        COMMAND ${Python3_EXECUTABLE} -I -B -S
            ${PROJECT_SOURCE_DIR}/bench/test_Wh2PublicPageAaR0.py)
    set_tests_properties(
        wirehair_wh2_public_page_aa_r0_selftest
        wirehair_wh2_public_page_aa_r0_python_selftest
        wirehair_wh2_public_page_aa_r0_python_test
        PROPERTIES TIMEOUT 120 WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        ENVIRONMENT "PYTHONWARNINGS=error"
        LABELS "wh2;benchmark;contract;public-api;page-aa")
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
       NOT CMAKE_CROSSCOMPILING AND BUILD_SHARED_LIBS)
        add_test(NAME wirehair_wh2_public_page_aa_r0_symbol_audit
            COMMAND ${CMAKE_COMMAND}
                "-DNM=${CMAKE_NM}"
                "-DPUBLIC_BORROWED_R1=$<TARGET_FILE:wirehair_wh2_public_page_aa_r0>"
                -P ${PROJECT_SOURCE_DIR}/cmake/Wh2PublicBorrowedCurrentVsWh1R1SymbolAudit.cmake)
        set_tests_properties(wirehair_wh2_public_page_aa_r0_symbol_audit
            PROPERTIES TIMEOUT 30 LABELS "wh2;benchmark;public-api;page-aa;symbols")
    endif()
    if(BUILD_CODEC_V2)
        # These native experiment targets execute test-hook-only equation
        # transforms.  The production library intentionally rejects them.
        add_executable(wirehair_wh2_native_codec_test
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h
            bench/Wh2NativeCodecTest.cpp)
        target_link_libraries(wirehair_wh2_native_codec_test PRIVATE
            wirehair_v2_policy)
        wirehair_apply_project_options(wirehair_wh2_native_codec_test)
        add_test(
            NAME wirehair_wh2_native_codec_test
            COMMAND wirehair_wh2_native_codec_test)
        set_tests_properties(wirehair_wh2_native_codec_test PROPERTIES
            TIMEOUT 180
            LABELS "wh2;benchmark;contract;native;codec")
        add_executable(wirehair_wh2_native_aa_calibration
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2NativeAaCalibration.cpp
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeAaCalibration.cpp")
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeAaCalibration.cpp"
            WIREHAIR_WH2_NATIVE_AA_SOURCE_SHA256)
        target_link_libraries(wirehair_wh2_native_aa_calibration PRIVATE
            wirehair_v2_timing_policy)
        target_compile_definitions(
            wirehair_wh2_native_aa_calibration PRIVATE
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\""
            "WIREHAIR_WH2_NATIVE_AA_SOURCE_SHA256=\"${WIREHAIR_WH2_NATIVE_AA_SOURCE_SHA256}\"")
        wirehair_apply_project_options(wirehair_wh2_native_aa_calibration)
        add_test(
            NAME wirehair_wh2_native_aa_calibration_selftest
            COMMAND wirehair_wh2_native_aa_calibration --selftest)
        set_tests_properties(
            wirehair_wh2_native_aa_calibration_selftest PROPERTIES
                TIMEOUT 180
                LABELS "wh2;benchmark;control;aa;native;timing")
        add_executable(wirehair_wh2_native_batch8_aa_calibration
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2NativeBatch8AaCalibration.cpp
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeBatch8AaCalibration.cpp")
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeBatch8AaCalibration.cpp"
            WIREHAIR_WH2_NATIVE_BATCH8_AA_SOURCE_SHA256)
        target_link_libraries(
            wirehair_wh2_native_batch8_aa_calibration PRIVATE
            wirehair_v2_timing_policy)
        target_compile_definitions(
            wirehair_wh2_native_batch8_aa_calibration PRIVATE
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\""
            "WIREHAIR_WH2_NATIVE_BATCH8_AA_SOURCE_SHA256=\"${WIREHAIR_WH2_NATIVE_BATCH8_AA_SOURCE_SHA256}\"")
        wirehair_apply_project_options(
            wirehair_wh2_native_batch8_aa_calibration)
        add_test(
            NAME wirehair_wh2_native_batch8_aa_calibration_selftest
            COMMAND wirehair_wh2_native_batch8_aa_calibration --selftest)
        set_tests_properties(
            wirehair_wh2_native_batch8_aa_calibration_selftest PROPERTIES
                TIMEOUT 180
                LABELS "wh2;benchmark;control;aa;batch8;native;timing")
        add_executable(wirehair_wh2_native_rdpru_batch128_aa_calibration
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2NativeRdpruAaCalibration.cpp
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeRdpruAaCalibration.cpp"
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeCodec.cpp"
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeCodec.h")
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeRdpruAaCalibration.cpp"
            WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256)
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeCodec.h"
            WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256)
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeCodec.cpp"
            WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256)
        file(SHA256
            "${PROJECT_SOURCE_DIR}/CMakeLists.txt"
            WIREHAIR_WH2_CMAKE_SHA256)
        target_link_libraries(
            wirehair_wh2_native_rdpru_batch128_aa_calibration PRIVATE
            wirehair_v2_timing_policy)
        target_compile_definitions(
            wirehair_wh2_native_rdpru_batch128_aa_calibration PRIVATE
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\""
            "WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256=\"${WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256}\""
            "WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256=\"${WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256}\""
            "WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256=\"${WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256}\""
            "WIREHAIR_WH2_CMAKE_SHA256=\"${WIREHAIR_WH2_CMAKE_SHA256}\"")
        wirehair_apply_project_options(
            wirehair_wh2_native_rdpru_batch128_aa_calibration)
        add_test(
            NAME wirehair_wh2_native_rdpru_batch128_aa_calibration_selftest
            COMMAND wirehair_wh2_native_rdpru_batch128_aa_calibration
                --selftest)
        set_tests_properties(
            wirehair_wh2_native_rdpru_batch128_aa_calibration_selftest
            PROPERTIES
                TIMEOUT 180
                LABELS "wh2;benchmark;control;aa;rdpru;batch128;native;counter")
        if(NOT "${WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256}" STREQUAL
           "0f553cce1df3acbc016250ccdf10fe3288980b3ae22bd4c6ee5c98e19bffa852")
            message(FATAL_ERROR
                "Frozen native RDPRU batch128 A/A v1 source hash changed")
        endif()
        add_executable(wirehair_wh2_native_rdpru_batch128_aa_calibration_v2
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2NativeRdpruAaCalibrationV2.cpp
            bench/Wh2RdpruTargetIdentityV2.cpp
            bench/Wh2RdpruTargetIdentityV2.h
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeRdpruAaCalibrationV2.cpp"
            "${PROJECT_SOURCE_DIR}/bench/Wh2RdpruTargetIdentityV2.cpp"
            "${PROJECT_SOURCE_DIR}/bench/Wh2RdpruTargetIdentityV2.h"
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeCodec.cpp"
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeCodec.h")
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2NativeRdpruAaCalibrationV2.cpp"
            WIREHAIR_WH2_NATIVE_RDPRU_AA_V2_SOURCE_SHA256)
        if(NOT
           "${WIREHAIR_WH2_NATIVE_RDPRU_AA_V2_SOURCE_SHA256}" STREQUAL
           "ac31b0c9b58e09b03a779c9c2d11205ff696ac3372b812b553822fe852be4ee8")
            message(FATAL_ERROR
                "Frozen native RDPRU batch128 A/A v2 source hash changed")
        endif()
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2RdpruTargetIdentityV2.cpp"
            WIREHAIR_WH2_RDPRU_TARGET_IDENTITY_V2_SOURCE_SHA256)
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2RdpruTargetIdentityV2.h"
            WIREHAIR_WH2_RDPRU_TARGET_IDENTITY_V2_HEADER_SHA256)
        if(NOT
           "${WIREHAIR_WH2_RDPRU_TARGET_IDENTITY_V2_SOURCE_SHA256}" STREQUAL
           "51ce3b59fcacd378e69da5c8e4a25327369a50fbae6809bccfeedcd4b72dc7b2")
            message(FATAL_ERROR
                "Frozen RDPRU target-identity v2 source hash changed")
        endif()
        if(NOT
           "${WIREHAIR_WH2_RDPRU_TARGET_IDENTITY_V2_HEADER_SHA256}" STREQUAL
           "af484c4382a5ec5ffab52fc3dfd5398b6d8dfed30d2347930c2f65c1787aa42c")
            message(FATAL_ERROR
                "Frozen RDPRU target-identity v2 header hash changed")
        endif()
        target_link_libraries(
            wirehair_wh2_native_rdpru_batch128_aa_calibration_v2 PRIVATE
            wirehair_v2_timing_policy)
        target_compile_definitions(
            wirehair_wh2_native_rdpru_batch128_aa_calibration_v2 PRIVATE
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\""
            "WIREHAIR_WH2_NATIVE_RDPRU_AA_V2_SOURCE_SHA256=\"${WIREHAIR_WH2_NATIVE_RDPRU_AA_V2_SOURCE_SHA256}\""
            "WIREHAIR_WH2_NATIVE_RDPRU_AA_V1_PARENT_SOURCE_SHA256=\"${WIREHAIR_WH2_NATIVE_RDPRU_AA_SOURCE_SHA256}\""
            "WIREHAIR_WH2_RDPRU_TARGET_IDENTITY_V2_SOURCE_SHA256=\"${WIREHAIR_WH2_RDPRU_TARGET_IDENTITY_V2_SOURCE_SHA256}\""
            "WIREHAIR_WH2_RDPRU_TARGET_IDENTITY_V2_HEADER_SHA256=\"${WIREHAIR_WH2_RDPRU_TARGET_IDENTITY_V2_HEADER_SHA256}\""
            "WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256=\"${WIREHAIR_WH2_NATIVE_CODEC_HEADER_SHA256}\""
            "WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256=\"${WIREHAIR_WH2_NATIVE_CODEC_SOURCE_SHA256}\""
            "WIREHAIR_WH2_CMAKE_SHA256=\"${WIREHAIR_WH2_CMAKE_SHA256}\"")
        wirehair_apply_project_options(
            wirehair_wh2_native_rdpru_batch128_aa_calibration_v2)
        add_test(
            NAME wirehair_wh2_native_rdpru_batch128_aa_calibration_v2_selftest
            COMMAND wirehair_wh2_native_rdpru_batch128_aa_calibration_v2
                --selftest)
        set_tests_properties(
            wirehair_wh2_native_rdpru_batch128_aa_calibration_v2_selftest
            PROPERTIES
                TIMEOUT 180
                LABELS
                    "wh2;benchmark;control;aa;rdpru;batch128;native;counter;target-v2")
        set(WIREHAIR_WH2_K3_HYPEROVAL_FALSIFIER_SOURCES
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2K3HyperovalCodec.cpp
            bench/Wh2K3HyperovalCodec.h
            bench/Wh2K3HyperovalFalsifier.cpp)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
            "${PROJECT_SOURCE_DIR}/bench/Wh2K3HyperovalCodec.cpp"
            "${PROJECT_SOURCE_DIR}/bench/Wh2K3HyperovalCodec.h"
            "${PROJECT_SOURCE_DIR}/bench/Wh2K3HyperovalFalsifier.cpp")
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2K3HyperovalCodec.cpp"
            WIREHAIR_WH2_K3_CODEC_SOURCE_SHA256)
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2K3HyperovalCodec.h"
            WIREHAIR_WH2_K3_CODEC_HEADER_SHA256)
        file(SHA256
            "${PROJECT_SOURCE_DIR}/bench/Wh2K3HyperovalFalsifier.cpp"
            WIREHAIR_WH2_K3_FALSIFIER_SOURCE_SHA256)
        set(WIREHAIR_WH2_K3_MANIFEST_DEFINITIONS
            "WIREHAIR_WH2_K3_LOW_CORPUS_SHA256=\"20cc19b3c3f46943cab8e9c7dcff74789b0417b2ff79d68c03cdab28acb4ce57\""
            "WIREHAIR_WH2_K3_HIGH_IID_SHA256=\"e03dab0bb1d0b891ac728aadde67b7e23b9bed0e1235bb0c6786bc4d46a65149\""
            "WIREHAIR_WH2_K3_SEQUENTIAL_SHA256=\"2124143403967032b04b1b707cacc819d9c812c08811a0bbd7546a1f49fe1b4c\""
            "WIREHAIR_WH2_K3_DELTA257_SHA256=\"990c36aaf9df5f8d4eb2d4e9dd5892b9389446af99e0e6c28de0fde190fbdd4d\""
            "WIREHAIR_WH2_K3_DELTA65535_SHA256=\"39e569a8bbabab6b013555cbf0d72c1352de08d6fea9d6ef226fad457b6f7ca7\""
            "WIREHAIR_WH2_K3_DESCENDING_SHA256=\"f86301a78dcfd6ca438cc0d77314cae46d2d15b2bb5da079f6c414b50616d2a1\""
            "WIREHAIR_WH2_K3_ONE_LOW_SHA256=\"de7e383e7e21f530105155388a9f26855c2658e0a239ec40837c969dbce311d0\""
            "WIREHAIR_WH2_K3_TWO_LOW_SHA256=\"4c95789ed24f3203f019ea462e6f40e5ea67d4fd20c6123e44cbdc6a7f964621\""
            "WIREHAIR_WH2_K3_CANONICAL_SHA256=\"d9b8b203d5cd21b3d846f5217fde0f7f0d3979b78992f51c30a96ebdee321d2e\""
            "WIREHAIR_WH2_K3_TRACE_SHA256=\"7cb9e6b5e6951eb84e65fa75faa279e927be1c455c6ce5f2cfd458f94aef9d93\""
            "WIREHAIR_WH2_K3_COMPLEMENT_SHA256=\"417aea72753f71a01b755fa414b59f655ea01922708a3001d9b0b455571981d9\"")
        add_executable(wirehair_wh2_k3_hyperoval_falsifier
            ${WIREHAIR_WH2_K3_HYPEROVAL_FALSIFIER_SOURCES})
        target_compile_definitions(
            wirehair_wh2_k3_hyperoval_falsifier PRIVATE
            ${WIREHAIR_WH2_K3_MANIFEST_DEFINITIONS}
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\""
            "WIREHAIR_WH2_K3_CODEC_SOURCE_SHA256=\"${WIREHAIR_WH2_K3_CODEC_SOURCE_SHA256}\""
            "WIREHAIR_WH2_K3_CODEC_HEADER_SHA256=\"${WIREHAIR_WH2_K3_CODEC_HEADER_SHA256}\""
            "WIREHAIR_WH2_K3_FALSIFIER_SOURCE_SHA256=\"${WIREHAIR_WH2_K3_FALSIFIER_SOURCE_SHA256}\"")
        target_link_libraries(
            wirehair_wh2_k3_hyperoval_falsifier PRIVATE
                ${WIREHAIR_INTERNAL_LIBRARY}
                Threads::Threads)
        wirehair_apply_project_options(
            wirehair_wh2_k3_hyperoval_falsifier)
        add_executable(wirehair_wh2_k3_hyperoval_hook_selftest
            ${WIREHAIR_WH2_K3_HYPEROVAL_FALSIFIER_SOURCES})
        target_compile_definitions(
            wirehair_wh2_k3_hyperoval_hook_selftest PRIVATE
            WIREHAIR_WH2_K3_HYPEROVAL_TEST_HOOKS=1
            ${WIREHAIR_WH2_K3_MANIFEST_DEFINITIONS}
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\""
            "WIREHAIR_WH2_K3_CODEC_SOURCE_SHA256=\"${WIREHAIR_WH2_K3_CODEC_SOURCE_SHA256}\""
            "WIREHAIR_WH2_K3_CODEC_HEADER_SHA256=\"${WIREHAIR_WH2_K3_CODEC_HEADER_SHA256}\""
            "WIREHAIR_WH2_K3_FALSIFIER_SOURCE_SHA256=\"${WIREHAIR_WH2_K3_FALSIFIER_SOURCE_SHA256}\"")
        target_link_libraries(
            wirehair_wh2_k3_hyperoval_hook_selftest PRIVATE
                ${WIREHAIR_INTERNAL_LIBRARY}
                Threads::Threads)
        wirehair_apply_project_options(
            wirehair_wh2_k3_hyperoval_hook_selftest)
        add_test(
            NAME wirehair_wh2_k3_hyperoval_falsifier_selftest
            COMMAND wirehair_wh2_k3_hyperoval_falsifier --selftest)
        add_test(
            NAME wirehair_wh2_k3_hyperoval_hook_selftest
            COMMAND wirehair_wh2_k3_hyperoval_hook_selftest --selftest)
        set_tests_properties(
            wirehair_wh2_k3_hyperoval_falsifier_selftest
            wirehair_wh2_k3_hyperoval_hook_selftest PROPERTIES
                TIMEOUT 120
                LABELS
                    "wh2;benchmark;contract;k3;hyperoval;falsifier")
        set(WIREHAIR_WH2_K2_PROJECTIVE_SCREEN_SOURCES
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2K2ProjectiveCodec.cpp
            bench/Wh2K2ProjectiveCodec.h
            bench/Wh2K2ProjectiveScreen.cpp
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h)
        add_executable(wirehair_wh2_k2_projective_screen
            ${WIREHAIR_WH2_K2_PROJECTIVE_SCREEN_SOURCES})
        target_link_libraries(wirehair_wh2_k2_projective_screen PRIVATE
            wirehair_v2_timing_policy)
        wirehair_apply_project_options(wirehair_wh2_k2_projective_screen)
        add_executable(wirehair_wh2_k2_projective_hook_selftest
            ${WIREHAIR_WH2_K2_PROJECTIVE_SCREEN_SOURCES})
        target_compile_definitions(
            wirehair_wh2_k2_projective_hook_selftest PRIVATE
            WIREHAIR_WH2_K2_PROJECTIVE_TEST_HOOKS=1)
        target_link_libraries(
            wirehair_wh2_k2_projective_hook_selftest PRIVATE
            wirehair_v2_timing_policy)
        wirehair_apply_project_options(
            wirehair_wh2_k2_projective_hook_selftest)
        add_test(
            NAME wirehair_wh2_k2_projective_screen_selftest
            COMMAND wirehair_wh2_k2_projective_screen --selftest)
        add_test(
            NAME wirehair_wh2_k2_projective_hook_selftest
            COMMAND wirehair_wh2_k2_projective_hook_selftest --selftest)
        set_tests_properties(
            wirehair_wh2_k2_projective_screen_selftest
            wirehair_wh2_k2_projective_hook_selftest PROPERTIES
                TIMEOUT 120
                LABELS
                    "wh2;benchmark;contract;native;k2;projective;falsifier")
        add_executable(wirehair_wh2_v2_facade_timing_worker_test
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2V2FacadeTimingWorker.cpp)
        target_link_libraries(
            wirehair_wh2_v2_facade_timing_worker_test PRIVATE wirehair)
        target_compile_definitions(
            wirehair_wh2_v2_facade_timing_worker_test PRIVATE
            WIREHAIR_WH2_FACADE_ROLE_CURRENT=1
            "WIREHAIR_WH2_HARNESS_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\""
            "WIREHAIR_WH2_IMPLEMENTATION_GIT_COMMIT=\"d6ddb35e046956174202584fb5a26c9a79679ea8\"")
        wirehair_apply_project_options(
            wirehair_wh2_v2_facade_timing_worker_test)
        add_test(
            NAME wirehair_wh2_v2_facade_timing_worker_selftest
            COMMAND wirehair_wh2_v2_facade_timing_worker_test --selftest)
        set_tests_properties(
            wirehair_wh2_v2_facade_timing_worker_selftest PROPERTIES
                TIMEOUT 120
                LABELS
                    "wh2;benchmark;contract;public-api;encoder;facade;timing")
        add_test(
            NAME wirehair_wh2_v2_facade_timing_screen_selftest
            COMMAND ${Python3_EXECUTABLE} -I -B
                ${PROJECT_SOURCE_DIR}/bench/Wh2V2FacadeTimingScreen.py
                --selftest)
        set_tests_properties(
            wirehair_wh2_v2_facade_timing_screen_selftest PROPERTIES
                TIMEOUT 60
                LABELS
                    "wh2;benchmark;contract;public-api;encoder;facade;timing")
        add_test(
            NAME wirehair_wh2_v2_facade_timing_launcher_selftest
            COMMAND ${Python3_EXECUTABLE} -I -B
                ${PROJECT_SOURCE_DIR}/bench/Wh2V2FacadeTimingLaunch.py
                --selftest)
        set_tests_properties(
            wirehair_wh2_v2_facade_timing_launcher_selftest PROPERTIES
                TIMEOUT 60
                LABELS
                    "wh2;benchmark;contract;public-api;encoder;facade;timing")
        add_test(
            NAME wirehair_wh2_v2_facade_timing_launcher_contract_test
            COMMAND ${Python3_EXECUTABLE} -I -B
                ${PROJECT_SOURCE_DIR}/bench/test_Wh2V2FacadeTimingLaunch.py)
        set_tests_properties(
            wirehair_wh2_v2_facade_timing_launcher_contract_test PROPERTIES
                TIMEOUT 120
                LABELS
                    "wh2;benchmark;contract;public-api;encoder;facade;timing")
        add_executable(wirehair_wh2_phase_attribution_test
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h
            bench/Wh2NativePanel.cpp
            bench/Wh2NativePanel.h
            bench/Wh2PhaseAttribution.cpp
            bench/Wh2PhaseAttribution.h
            bench/Wh2PhaseAttributionTest.cpp)
        target_link_libraries(wirehair_wh2_phase_attribution_test PRIVATE
            wirehair_v2_timing_policy)
        wirehair_apply_project_options(wirehair_wh2_phase_attribution_test)
        add_test(
            NAME wirehair_wh2_phase_attribution_test
            COMMAND wirehair_wh2_phase_attribution_test)
        set_tests_properties(wirehair_wh2_phase_attribution_test PROPERTIES
            TIMEOUT 180
            LABELS "wh2;benchmark;contract;native;timing;phase")
        add_executable(wirehair_wh2_systematic_emission_screen
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h
            bench/Wh2NativePanel.cpp
            bench/Wh2NativePanel.h
            bench/Wh2SystematicEmissionScreen.cpp)
        target_link_libraries(wirehair_wh2_systematic_emission_screen PRIVATE
            wirehair_v2_timing_policy)
        target_compile_definitions(
            wirehair_wh2_systematic_emission_screen PRIVATE
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\"")
        wirehair_apply_project_options(
            wirehair_wh2_systematic_emission_screen)
        add_test(
            NAME wirehair_wh2_systematic_emission_screen_selftest
            COMMAND wirehair_wh2_systematic_emission_screen --selftest)
        set_tests_properties(
            wirehair_wh2_systematic_emission_screen_selftest PROPERTIES
                TIMEOUT 30
                LABELS "wh2;benchmark;contract;native;encoder;systematic")
        add_test(
            NAME wirehair_wh2_systematic_emission_analyzer_selftest
            COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/bench/Wh2SystematicEmissionScreen.py
                --selftest)
        set_tests_properties(
            wirehair_wh2_systematic_emission_analyzer_selftest PROPERTIES
                TIMEOUT 30
                LABELS "wh2;benchmark;contract;native;encoder;systematic")
        add_executable(wirehair_wh2_broad_systematic_emission_screen
            bench/Wh2BroadSystematicEmissionScreen.cpp
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h
            bench/Wh2NativePanel.cpp
            bench/Wh2NativePanel.h)
        target_link_libraries(
            wirehair_wh2_broad_systematic_emission_screen PRIVATE
            wirehair_v2_timing_policy)
        target_compile_definitions(
            wirehair_wh2_broad_systematic_emission_screen PRIVATE
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\"")
        wirehair_apply_project_options(
            wirehair_wh2_broad_systematic_emission_screen)
        add_test(
            NAME wirehair_wh2_broad_systematic_emission_screen_selftest
            COMMAND wirehair_wh2_broad_systematic_emission_screen --selftest)
        set_tests_properties(
            wirehair_wh2_broad_systematic_emission_screen_selftest PROPERTIES
                TIMEOUT 120
                LABELS
                    "wh2;benchmark;contract;native;encoder;systematic;broad")
        add_test(
            NAME wirehair_wh2_broad_systematic_emission_analyzer_selftest
            COMMAND ${Python3_EXECUTABLE} -I -B
                ${PROJECT_SOURCE_DIR}/bench/Wh2BroadSystematicEmissionScreen.py
                --selftest)
        set_tests_properties(
            wirehair_wh2_broad_systematic_emission_analyzer_selftest PROPERTIES
                TIMEOUT 30
                LABELS
                    "wh2;benchmark;contract;native;encoder;systematic;broad")
        if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
           NOT CMAKE_CROSSCOMPILING AND NOT CMAKE_NM)
            message(FATAL_ERROR
                "native Linux broad systematic-emission contract requires CMAKE_NM audits")
        endif()
        if(CMAKE_NM)
            add_test(
                NAME wirehair_wh2_broad_systematic_emission_symbol_audit
                COMMAND ${CMAKE_COMMAND}
                    "-DNM=${CMAKE_NM}"
                    "-DBROAD_SYSTEMATIC_EMISSION=$<TARGET_FILE:wirehair_wh2_broad_systematic_emission_screen>"
                    -P ${PROJECT_SOURCE_DIR}/cmake/Wh2BroadSystematicEmissionSymbolAudit.cmake)
            set_tests_properties(
                wirehair_wh2_broad_systematic_emission_symbol_audit PROPERTIES
                    TIMEOUT 30
                    LABELS
                        "wh2;benchmark;contract;native;encoder;systematic;broad;symbols")
        endif()
        add_executable(wirehair_wh2_direct_systematic_complement_screen
            bench/Wh2DirectSystematicComplementScreen.cpp
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h
            bench/Wh2NativePanel.cpp
            bench/Wh2NativePanel.h)
        target_link_libraries(
            wirehair_wh2_direct_systematic_complement_screen PRIVATE
            wirehair_v2_timing_policy)
        target_compile_definitions(
            wirehair_wh2_direct_systematic_complement_screen PRIVATE
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\"")
        wirehair_apply_project_options(
            wirehair_wh2_direct_systematic_complement_screen)
        add_test(
            NAME wirehair_wh2_direct_systematic_complement_screen_selftest
            COMMAND wirehair_wh2_direct_systematic_complement_screen
                --selftest)
        set_tests_properties(
            wirehair_wh2_direct_systematic_complement_screen_selftest
            PROPERTIES
                TIMEOUT 120
                LABELS
                    "wh2;benchmark;contract;native;encoder;systematic;complement")
        if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
           NOT CMAKE_CROSSCOMPILING)
            add_test(
                NAME wirehair_wh2_direct_systematic_complement_analyzer_selftest
                COMMAND ${Python3_EXECUTABLE} -I -B
                    ${PROJECT_SOURCE_DIR}/bench/Wh2DirectSystematicComplementScreen.py
                    --selftest)
            set_tests_properties(
                wirehair_wh2_direct_systematic_complement_analyzer_selftest
                PROPERTIES
                    TIMEOUT 30
                    LABELS
                        "wh2;benchmark;contract;native;encoder;systematic;complement")
            add_test(
                NAME wirehair_wh2_direct_systematic_complement_launcher_selftest
                COMMAND ${Python3_EXECUTABLE} -I -B
                    ${PROJECT_SOURCE_DIR}/bench/test_wh2_direct_systematic_complement_launch.py)
            set_tests_properties(
                wirehair_wh2_direct_systematic_complement_launcher_selftest
                PROPERTIES
                    TIMEOUT 30
                    LABELS
                        "wh2;benchmark;contract;native;encoder;systematic;complement;launcher")
            add_test(
                NAME wirehair_expo_thermal_sampler_test
                COMMAND ${CMAKE_COMMAND} -E env
                    PYTHON=${Python3_EXECUTABLE}
                    bash
                    ${PROJECT_SOURCE_DIR}/bench/test_wirehair_expo_thermal_sampler.sh)
            set_tests_properties(wirehair_expo_thermal_sampler_test
                PROPERTIES
                    LABELS
                        "wh2;benchmark;contract;complement;v2;research;thermal;plausibility;transactional"
                    TIMEOUT 30)
        endif()
        if(CMAKE_NM)
            add_test(
                NAME wirehair_wh2_direct_systematic_complement_symbol_audit
                COMMAND ${CMAKE_COMMAND}
                    "-DNM=${CMAKE_NM}"
                    "-DDIRECT_SYSTEMATIC_COMPLEMENT=$<TARGET_FILE:wirehair_wh2_direct_systematic_complement_screen>"
                    -P ${PROJECT_SOURCE_DIR}/cmake/Wh2DirectSystematicComplementSymbolAudit.cmake)
            set_tests_properties(
                wirehair_wh2_direct_systematic_complement_symbol_audit
                PROPERTIES
                    TIMEOUT 30
                    LABELS
                        "wh2;benchmark;contract;native;encoder;systematic;complement;symbols")
        endif()
        if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
           NOT CMAKE_CROSSCOMPILING)
            add_executable(wirehair_wh2_mix2_seed_repair_worker
                bench/Wh2Mix2SeedRepairWorker.cpp
                bench/Wh2FrozenTrace.cpp
                bench/Wh2FrozenTrace.h
                bench/Wh2NativeCodec.cpp
                bench/Wh2NativeCodec.h)
            target_compile_definitions(
                wirehair_wh2_mix2_seed_repair_worker PRIVATE
                "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\"")
            target_link_libraries(
                wirehair_wh2_mix2_seed_repair_worker PRIVATE
                wirehair_v2_timing_policy
                Threads::Threads)
            wirehair_apply_project_options(
                wirehair_wh2_mix2_seed_repair_worker)
            add_test(
                NAME wirehair_wh2_mix2_seed_repair_worker_selftest
                COMMAND wirehair_wh2_mix2_seed_repair_worker --self-test)
            set_tests_properties(
                wirehair_wh2_mix2_seed_repair_worker_selftest PROPERTIES
                    TIMEOUT 120
                    LABELS
                        "wh2;benchmark;contract;native;recovery;mix2;repair")
            add_test(
                NAME wirehair_wh2_mix2_seed_repair_worker_controller_test
                COMMAND ${Python3_EXECUTABLE} -B
                    ${PROJECT_SOURCE_DIR}/bench/test_wh2_mix2_seed_repair_worker.py
                    $<TARGET_FILE:wirehair_wh2_mix2_seed_repair_worker>)
            set_tests_properties(
                wirehair_wh2_mix2_seed_repair_worker_controller_test
                PROPERTIES
                    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                    ENVIRONMENT "PYTHONWARNINGS=error"
                    TIMEOUT 120
                    LABELS
                        "wh2;benchmark;contract;native;recovery;mix2;repair")
            add_test(
                NAME wirehair_wh2_mix2_promotion_short_screen_test
                COMMAND ${Python3_EXECUTABLE} -B
                    ${PROJECT_SOURCE_DIR}/bench/test_wh2_mix2_promotion_short_screen.py)
            set_tests_properties(
                wirehair_wh2_mix2_promotion_short_screen_test PROPERTIES
                    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                    ENVIRONMENT "PYTHONWARNINGS=error"
                    TIMEOUT 30
                    LABELS
                        "wh2;benchmark;contract;native;recovery;mix2;promotion")
        endif()
        add_executable(wirehair_wh2_contract_worker
            bench/Wh2ContractWorker.cpp
            bench/Wh2FrozenTrace.cpp
            bench/Wh2FrozenTrace.h
            bench/Wh2FrozenTiming.cpp
            bench/Wh2FrozenTiming.h
            bench/Wh2NativeCodec.cpp
            bench/Wh2NativeCodec.h
            bench/Wh2NativePanel.cpp
            bench/Wh2NativePanel.h
            bench/Wh2PhaseAttribution.cpp
            bench/Wh2PhaseAttribution.h)
        target_compile_definitions(wirehair_wh2_contract_worker PRIVATE
            "WIREHAIR_WH2_SOURCE_GIT_COMMIT=\"${WIREHAIR_WH2_SOURCE_GIT_COMMIT}\"")
        target_link_libraries(wirehair_wh2_contract_worker PRIVATE
            wirehair_v2_timing_policy
            Threads::Threads)
        wirehair_apply_project_options(wirehair_wh2_contract_worker)
        if(CMAKE_NM)
            add_test(
                NAME wirehair_wh2_timing_policy_symbol_audit
                COMMAND ${CMAKE_COMMAND}
                    "-DNM=${CMAKE_NM}"
                    "-DOBJDUMP=${CMAKE_OBJDUMP}"
                    "-DSYSTEM_NAME=${CMAKE_SYSTEM_NAME}"
                    "-DLTO_MODE=${WH_LTO}"
                    "-DTIMING_POLICY=$<TARGET_FILE:wirehair_v2_timing_policy>"
                    "-DCONTRACT_WORKER=$<TARGET_FILE:wirehair_wh2_contract_worker>"
                    "-DPHASE_TIMING=$<TARGET_FILE:wirehair_wh2_phase_attribution_test>"
                    "-DSYSTEMATIC_EMISSION=$<TARGET_FILE:wirehair_wh2_systematic_emission_screen>"
                    "-DBROAD_SYSTEMATIC_EMISSION=$<TARGET_FILE:wirehair_wh2_broad_systematic_emission_screen>"
                    "-DDIRECT_SYSTEMATIC_COMPLEMENT=$<TARGET_FILE:wirehair_wh2_direct_systematic_complement_screen>"
                    -P ${PROJECT_SOURCE_DIR}/cmake/Wh2TimingPolicySymbolAudit.cmake)
            set_tests_properties(
                wirehair_wh2_timing_policy_symbol_audit
                PROPERTIES
                    TIMEOUT 30
                    LABELS "wh2;benchmark;contract;timing;symbols")
        endif()
    endif()
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND NOT CMAKE_CROSSCOMPILING)
        if(BUILD_CODEC_V2)
            add_test(
                NAME wirehair_wh2_contract_worker_cli_test
                COMMAND ${Python3_EXECUTABLE}
                    ${PROJECT_SOURCE_DIR}/bench/wh2_contract_worker_cli_test.py
                    $<TARGET_FILE:wirehair_wh2_contract_worker>)
            set_tests_properties(wirehair_wh2_contract_worker_cli_test
                PROPERTIES
                WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                ENVIRONMENT "PYTHONWARNINGS=error"
                TIMEOUT 120
                LABELS "wh2;benchmark;contract;native;worker;recovery")
        endif()
        add_test(
            NAME wh2_run_native_recovery_screen_test
            COMMAND ${Python3_EXECUTABLE} -m unittest -v
                bench/test_wh2_run_native_recovery_screen.py)
        set_tests_properties(wh2_run_native_recovery_screen_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 180
            LABELS "wh2;benchmark;contract;native;runner;recovery")
        add_test(
            NAME wh2_precodefail_work_screen_test
            COMMAND ${Python3_EXECUTABLE} -m unittest -v
                bench/test_wh2_precodefail_work_screen.py)
        set_tests_properties(wh2_precodefail_work_screen_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 180
            LABELS "wh2;benchmark;contract;recovery;precodefail;runner")
        add_test(
            NAME wh2_mix2_rank_work_screen_test
            COMMAND ${Python3_EXECUTABLE} -B -m unittest -v
                bench/test_wh2_mix2_rank_work_screen.py)
        set_tests_properties(wh2_mix2_rank_work_screen_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS "wh2;benchmark;contract;recovery;precodefail;runner")
        add_test(
            NAME wh2_mix2_attempt_crossfit_test
            COMMAND ${Python3_EXECUTABLE} -B -m unittest -v
                bench/test_wh2_mix2_attempt_crossfit.py)
        set_tests_properties(wh2_mix2_attempt_crossfit_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS
                "wh2;benchmark;contract;recovery;precodefail;mix2;selection")
        add_test(
            NAME wh2_mix2_pair_stage0_test
            COMMAND ${Python3_EXECUTABLE} -B -m unittest -v
                bench/test_wh2_mix2_pair_stage0.py)
        set_tests_properties(wh2_mix2_pair_stage0_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS
                "wh2;benchmark;contract;recovery;precodefail;mix2;selection")
        add_test(
            NAME wh2_mix2_packet_offset_stage0_test
            COMMAND ${Python3_EXECUTABLE} -B -m unittest -v
                bench/test_wh2_mix2_packet_offset_stage0.py)
        set_tests_properties(wh2_mix2_packet_offset_stage0_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS
                "wh2;benchmark;contract;recovery;precodefail;mix2;selection")
        add_test(
            NAME wh2_mix2_pair_stage1_test
            COMMAND ${Python3_EXECUTABLE} -B -m unittest -v
                bench/test_wh2_mix2_pair_stage1.py)
        set_tests_properties(wh2_mix2_pair_stage1_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS
                "wh2;benchmark;contract;recovery;precodefail;mix2;selection")
        add_test(
            NAME wh2_mix2_packet_offset_stage1_test
            COMMAND ${Python3_EXECUTABLE} -B -m unittest -v
                bench/test_wh2_mix2_packet_offset_stage1.py)
        set_tests_properties(wh2_mix2_packet_offset_stage1_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS
                "wh2;benchmark;contract;recovery;precodefail;mix2;selection")
        add_test(
            NAME wh2_mix2_packet_offset_stage2_test
            COMMAND ${Python3_EXECUTABLE} -B -m unittest -v
                bench/test_wh2_mix2_packet_offset_stage2.py)
        set_tests_properties(wh2_mix2_packet_offset_stage2_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS
                "wh2;benchmark;contract;recovery;precodefail;mix2;selection")
        add_test(
            NAME wh2_mix2_packet_offset_stage3_test
            COMMAND ${Python3_EXECUTABLE} -B -m unittest -v
                bench/test_wh2_mix2_packet_offset_stage3.py)
        set_tests_properties(wh2_mix2_packet_offset_stage3_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS
                "wh2;benchmark;contract;recovery;precodefail;mix2;selection")
        add_test(
            NAME wh2_mix2_seed_repair_contract_test
            COMMAND ${Python3_EXECUTABLE} -B -m unittest -v
                bench/test_wh2_mix2_seed_repair_contract.py)
        set_tests_properties(wh2_mix2_seed_repair_contract_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS "wh2;benchmark;contract;native;recovery;mix2;repair")
        add_test(
            NAME wh2_select_development_architecture_test
            COMMAND ${Python3_EXECUTABLE} -m unittest -v
                bench/test_wh2_select_development_architecture.py)
        set_tests_properties(wh2_select_development_architecture_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 120
            LABELS "wh2;benchmark;contract;selection;recovery;timing")
        add_test(
            NAME wh2_phase_attribution_python_test
            COMMAND ${Python3_EXECUTABLE} -m unittest -v
                bench/test_wh2_phase_attribution.py)
        set_tests_properties(wh2_phase_attribution_python_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 180
            LABELS "wh2;benchmark;contract;native;timing;phase")
        add_test(
            NAME wh2_run_native_phase_attribution_test
            COMMAND ${Python3_EXECUTABLE} -m unittest -v
                bench/test_wh2_run_native_phase_attribution.py)
        set_tests_properties(wh2_run_native_phase_attribution_test PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            ENVIRONMENT "PYTHONWARNINGS=error"
            TIMEOUT 180
            LABELS "wh2;benchmark;contract;native;runner;timing;phase")
    endif()
    add_test(
        NAME wh2_benchmark_contract_test
        COMMAND ${Python3_EXECUTABLE} -m unittest -v
            bench/test_wh2_benchmark_contract.py)
    set_tests_properties(wh2_benchmark_contract_test PROPERTIES
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        ENVIRONMENT "PYTHONWARNINGS=error"
        TIMEOUT 180
        LABELS "wh2;benchmark;contract")
    add_test(
        NAME wh2_native_short_screen_test
        COMMAND ${Python3_EXECUTABLE} -m unittest -v
            bench/test_wh2_native_short_screen.py)
    set_tests_properties(wh2_native_short_screen_test PROPERTIES
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        ENVIRONMENT "PYTHONWARNINGS=error"
        TIMEOUT 360
        LABELS "wh2;benchmark;contract;native")
    add_test(
        NAME wh2_run_native_short_screen_test
        COMMAND ${Python3_EXECUTABLE} -m unittest -v
            bench/test_wh2_run_native_short_screen.py)
    set_tests_properties(wh2_run_native_short_screen_test PROPERTIES
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        ENVIRONMENT "PYTHONWARNINGS=error"
        # The fault-injection/cleanup matrix is validated separately on Python
        # 3.8 and current interpreters and takes roughly five to six minutes
        # per configured interpreter on one core under full-system load.
        TIMEOUT 600
        LABELS "wh2;benchmark;contract;native;runner;affinity")
endif()

if(BUILD_TESTS AND NOT CMAKE_CROSSCOMPILING)
    find_package(PkgConfig QUIET)
    get_property(WIREHAIR_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    foreach(_wirehair_package_variant IN ITEMS static shared both)
        if(_wirehair_package_variant STREQUAL "shared")
            set(_wirehair_package_shared ON)
            set(_wirehair_package_both OFF)
            set(_wirehair_package_config Release)
        elseif(_wirehair_package_variant STREQUAL "both")
            set(_wirehair_package_shared OFF)
            set(_wirehair_package_both ON)
            set(_wirehair_package_config RelWithDebInfo)
        else()
            set(_wirehair_package_shared OFF)
            set(_wirehair_package_both OFF)
            set(_wirehair_package_config Debug)
        endif()
        add_test(
            NAME package_${_wirehair_package_variant}_e2e
            COMMAND ${CMAKE_COMMAND}
                -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/package-${_wirehair_package_variant}
                -DTEST_GENERATOR=${CMAKE_GENERATOR}
                -DTEST_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
                -DTEST_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}
                -DTEST_C_COMPILER=${CMAKE_C_COMPILER}
                -DTEST_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DTEST_SHARED=${_wirehair_package_shared}
                -DTEST_BOTH=${_wirehair_package_both}
                -DTEST_CONFIG=${_wirehair_package_config}
                -DTEST_STRICT=${WIREHAIR_STRICT_WARNINGS}
                -DTEST_NM=${CMAKE_NM}
                -DTEST_PKG_CONFIG=${PKG_CONFIG_EXECUTABLE}
                -DTEST_PROJECT_VERSION=${PROJECT_VERSION}
                -DTEST_PYTHON_EXECUTABLE=${Python3_EXECUTABLE}
                -P ${PROJECT_SOURCE_DIR}/test/package/RunPackageE2E.cmake)
        set_tests_properties(package_${_wirehair_package_variant}_e2e
            PROPERTIES
                TIMEOUT 240
                LABELS "package;cmake;e2e;ci-once")
    endforeach()
    if(UNIX AND NOT MSVC)
        if(APPLE)
            set(WIREHAIR_LINKER_SENTINEL "-Wl,-undefined,error")
            set(WIREHAIR_LINKER_SENTINEL_PATTERN "-undefined")
            set(WIREHAIR_PIC_FLAG_REQUIRED OFF)
        else()
            set(WIREHAIR_LINKER_SENTINEL "-Wl,--no-undefined")
            set(WIREHAIR_LINKER_SENTINEL_PATTERN "--no-undefined")
            set(WIREHAIR_PIC_FLAG_REQUIRED ON)
        endif()
        add_test(
            NAME build_policy_e2e
            COMMAND ${CMAKE_COMMAND}
                -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/build-policy
                -DTEST_GENERATOR=${CMAKE_GENERATOR}
                -DTEST_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
                -DTEST_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}
                -DTEST_C_COMPILER=${CMAKE_C_COMPILER}
                -DTEST_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DTEST_CONFIG=Release
                -DTEST_MULTI_CONFIG=${WIREHAIR_MULTI_CONFIG}
                -DTEST_EXE_SUFFIX=${CMAKE_EXECUTABLE_SUFFIX}
                -DTEST_LINKER_SENTINEL=${WIREHAIR_LINKER_SENTINEL}
                -DTEST_LINKER_PATTERN=${WIREHAIR_LINKER_SENTINEL_PATTERN}
                -DTEST_NATIVE_SUPPORTED=${WIREHAIR_HAS_MARCH_NATIVE}
                -DTEST_PIC_FLAG_REQUIRED=${WIREHAIR_PIC_FLAG_REQUIRED}
                -P ${PROJECT_SOURCE_DIR}/test/cmake/RunBuildPolicyE2E.cmake)
        set_tests_properties(build_policy_e2e PROPERTIES
            TIMEOUT 240
            RUN_SERIAL TRUE
            LABELS "build;cmake;e2e;ci-once")
    endif()
endif()
