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
option(BUILD_CODEC_V2 "Build the isolated v2 codec scaffold" ON)
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
option(WIREHAIR_STRICT_WARNINGS
    "Treat compiler warnings as errors for Wirehair targets" OFF)
set(WH_LTO "OFF" CACHE STRING "Link-time optimization mode: OFF, ON, AUTO, or THIN")
set_property(CACHE WH_LTO PROPERTY STRINGS OFF ON AUTO THIN)
set(WH_PGO_MODE "OFF" CACHE STRING "Profile-guided optimization mode: OFF, GENERATE, or USE")
set_property(CACHE WH_PGO_MODE PROPERTY STRINGS OFF GENERATE USE)
set(WH_PGO_DIR "${CMAKE_BINARY_DIR}/pgo" CACHE PATH "Profile data directory for WH_PGO_MODE")

set(LIB_SOURCE_FILES
        wirehair.cpp
        include/wirehair/wirehair.h
        gf256.cpp
        gf256.h
        WirehairCodec.cpp
        WirehairCodec.h
        WirehairCompiler.h
        WirehairEnvironment.h
        WirehairTools.cpp
        WirehairTools.h
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
    target_link_libraries(${target} PRIVATE Threads::Threads)
    if(variant STREQUAL "SHARED")
        target_compile_definitions(${target}
            PRIVATE WIREHAIR_BUILDING=1 WIREHAIR_DLL=1
            INTERFACE WIREHAIR_DLL=1)
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
                OUTPUT_NAME wirehair_shared)
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
if(WIN32 AND WIREHAIR_PRIMARY_VARIANT STREQUAL "SHARED")
    if(TARGET wirehair_static)
        set(WIREHAIR_INTERNAL_LIBRARY wirehair_static)
    else()
        # The DLL exposes the public C ABI only; V2 and offline generators use
        # internal C++ codec types.  Their support stays out of routine builds
        # unless an internal target needs it.
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
install(FILES python/whirehair.py DESTINATION python)

set(WIREHAIR_PACKAGE_INSTALL_DIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/wirehair")
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

if(BUILD_TESTS AND NOT CMAKE_CROSSCOMPILING)
    get_property(WIREHAIR_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    foreach(_wirehair_package_variant IN ITEMS static shared)
        if(_wirehair_package_variant STREQUAL "shared")
            set(_wirehair_package_shared ON)
            set(_wirehair_package_config Release)
        else()
            set(_wirehair_package_shared OFF)
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
                -DTEST_CONFIG=${_wirehair_package_config}
                -DTEST_STRICT=${WIREHAIR_STRICT_WARNINGS}
                -P ${PROJECT_SOURCE_DIR}/test/package/RunPackageE2E.cmake)
        set_tests_properties(package_${_wirehair_package_variant}_e2e
            PROPERTIES TIMEOUT 240)
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
            LABELS "build;cmake;e2e")
    endif()
endif()
