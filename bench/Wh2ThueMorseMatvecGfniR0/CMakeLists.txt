cmake_minimum_required(VERSION 3.16)
project(Wh2ThueMorseMatvecGfniR0 LANGUAGES CXX)

# Standalone neutral arithmetic tests only; no codec, lookup or timing target.
option(WH2_MATVEC_SANITIZERS "Enable ASan/UBSan for neutral tests" OFF)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_CROSSCOMPILING OR
   NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CONFIGURATION_TYPES)
    message(FATAL_ERROR "Native single-configuration Linux/GNU build required")
endif()
get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." REALPATH)
get_filename_component(_binary "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)
string(FIND "${_binary}/" "${_root}/" _inside)
if(_inside EQUAL 0)
    message(FATAL_ERROR "Use a fresh external build directory")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" _config)
set(_flags "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_config}} ${CMAKE_EXE_LINKER_FLAGS}")
if(CMAKE_INTERPROCEDURAL_OPTIMIZATION OR CMAKE_INTERPROCEDURAL_OPTIMIZATION_${_config} OR
   _flags MATCHES "WH_COUNT|WIREHAIR_TESTING|WIREHAIR_V2_ENABLE_TEST_HOOKS|-DGF256|-flto|-fprofile|-m[a-z]")
    message(FATAL_ERROR "Counters, hooks, LTO, PGO and injected ISA flags are forbidden")
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(Threads REQUIRED)
enable_testing()

# The original object is built ONLY for a machine-code comparison. Linking it
# beside the candidate would duplicate the shared GF context and is forbidden.
add_library(matvec_original_gf_audit OBJECT "${_root}/gf256.cpp")
add_library(matvec_candidate OBJECT "${_root}/bench/Wh2ThueMorseMatvecGfniR0.cpp")
add_executable(wh2_matvec_gfni_test
    "${_root}/bench/Wh2ThueMorseMatvecGfniR0Test.cpp"
    $<TARGET_OBJECTS:matvec_candidate>)
foreach(_target IN ITEMS matvec_original_gf_audit matvec_candidate wh2_matvec_gfni_test)
    target_include_directories(${_target} PRIVATE "${_root}" "${_root}/bench")
    set_target_properties(${_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_compile_options(${_target} PRIVATE -Wall -Wextra -Wpedantic -Werror -fno-lto)
    if(WH2_MATVEC_SANITIZERS)
        target_compile_options(${_target} PRIVATE -O1 -g -fno-omit-frame-pointer
            -fsanitize=address,undefined -fno-pie)
    else()
        target_compile_options(${_target} PRIVATE -O3)
    endif()
endforeach()
target_link_libraries(wh2_matvec_gfni_test PRIVATE Threads::Threads)
target_link_options(wh2_matvec_gfni_test PRIVATE -fno-lto
    "-Wl,-Map,${_binary}/matvec-test.map")
if(WH2_MATVEC_SANITIZERS)
    target_link_options(wh2_matvec_gfni_test PRIVATE -fsanitize=address,undefined -no-pie)
endif()
add_test(NAME matvec_native COMMAND wh2_matvec_gfni_test --require-gfni)
set_tests_properties(matvec_native PROPERTIES TIMEOUT 60)
if(WH2_MATVEC_SANITIZERS)
    # Exercise the fake-stack layout that exposed an aligned initializer store.
    set_tests_properties(matvec_native PROPERTIES ENVIRONMENT
        "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:detect_stack_use_after_return=1;UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1")
endif()
if(NOT WH2_MATVEC_SANITIZERS)
    find_program(_qemu qemu-x86_64)
    if(_qemu)
        add_test(NAME matvec_haswell COMMAND "${_qemu}" -cpu Haswell
            $<TARGET_FILE:wh2_matvec_gfni_test> --expect-scalar)
        set_tests_properties(matvec_haswell PROPERTIES TIMEOUT 60)
    endif()
endif()
