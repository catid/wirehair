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
option(WIREHAIR_ENABLE_LIBFUZZER
    "Build V2 fuzz entry points with Clang libFuzzer main/coverage" OFF)
option(WIREHAIR_STRICT_WARNINGS
    "Treat compiler warnings as errors for Wirehair targets" OFF)
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
        codec/WirehairV2GF16.cpp
        codec/WirehairV2GF16.h
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
        codec/WirehairV2MixedBuckets.h
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

if(BUILD_TESTS AND NOT CMAKE_CROSSCOMPILING)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    find_package(PkgConfig QUIET)
    find_program(WIREHAIR_GIT_EXECUTABLE git)
    # The campaign stack imports fcntl and these suites exercise POSIX
    # path/executable semantics.  Keep them off native Windows rather than
    # registering tests that fail during module import.
    if(Python3_Interpreter_FOUND AND UNIX)
        foreach(_wh2_preferred_suite IN ITEMS campaign holdout timing)
            add_test(
                NAME wh2_preferred_attempt_${_wh2_preferred_suite}_test
                COMMAND ${Python3_EXECUTABLE}
                    ${PROJECT_SOURCE_DIR}/bench/wh2_preferred_attempt_${_wh2_preferred_suite}_test.py
                    -v)
            set_tests_properties(
                wh2_preferred_attempt_${_wh2_preferred_suite}_test
                PROPERTIES
                    ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/bench"
                    LABELS
                        "v2;research;state-machine;${_wh2_preferred_suite}"
                    TIMEOUT 120)
        endforeach()
        add_test(
            NAME wh2_rank_floor_two_anchor_screen_test
            COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/bench/wh2_rank_floor_two_anchor_screen_test.py
                -v)
        set_tests_properties(wh2_rank_floor_two_anchor_screen_test
            PROPERTIES
                ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/bench"
                LABELS "v2;research;thermal"
                TIMEOUT 30)
    endif()
    if(Python3_Interpreter_FOUND AND WIREHAIR_GIT_EXECUTABLE AND UNIX)
        add_test(
            NAME wh2_preferred_attempt_search_test
            COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/bench/wh2_preferred_attempt_search_test.py
                -v)
        set_tests_properties(wh2_preferred_attempt_search_test
            PROPERTIES
                ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/bench"
                LABELS "v2;research;state-machine;beacon"
                TIMEOUT 60)
    endif()
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
