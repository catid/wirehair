cmake_minimum_required(VERSION 3.20)
project(Wh2K6Serialized LANGUAGES C CXX)
option(WH2_K6_SANITIZERS "ASan/UBSan neutral qualification" OFF)
option(WH2_K6_SCALAR "Portable GF backend neutral qualification on this host" OFF)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_CROSSCOMPILING OR
   NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CONFIGURATION_TYPES)
    message(FATAL_ERROR "This isolated qualification build requires native Linux/GNU")
endif()
get_filename_component(_root "${CMAKE_CURRENT_SOURCE_DIR}/../.." REALPATH)
get_filename_component(_binary "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)
string(FIND "${_binary}/" "${_root}/" _inside)
if(_inside EQUAL 0)
    message(FATAL_ERROR "Use an external build directory")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" _config)
set(_flags "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${_config}} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_config}} ${CMAKE_EXE_LINKER_FLAGS}")
if(CMAKE_INTERPROCEDURAL_OPTIMIZATION OR CMAKE_INTERPROCEDURAL_OPTIMIZATION_${_config} OR
   _flags MATCHES "WH_COUNT|WIREHAIR_TESTING|WIREHAIR_V2_ENABLE_TEST_HOOKS|-DGF256|-flto|-fprofile|-m[a-z]|-fshort-enums")
    message(FATAL_ERROR "Counters, hooks, LTO, PGO, injected ISA and enum ABI flags are forbidden")
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(Threads REQUIRED)
enable_testing()

set(_pins
    "bench/Wh2ThueMorseNativeCodec.cpp|8e8e32e05394c16af0d599d83db4d319aa429d6c6f0377714cdfda540ff60c2d"
    "bench/Wh2ThueMorseNativeCodec.h|0fedb1d611d75a1235abf92b99761c5daee7ecb5b43efe8007c0f5a5113b037b"
    "bench/Wh2ThueMorseTinyGfniR0Codec.cpp|0609762920b7ddfd108a77473da2538714aecd8e65427bc624a4d3659cf8d71e"
    "bench/Wh2ThueMorseTinyGfniR0.cpp|d596555cefc7486f4a101f9bcd10b4ac43f634bdbdb91e3d0f85619f229b54ea"
    "bench/Wh2ThueMorseTinyGfniR0.h|3c562e03ea7e90061a80c598c22cf8060fa2f32151c7a2fd7120d8764df39068"
    "bench/Wh2ThueMorseTinyPayloadR0.h|e38c818bfe3be9aa4a7f4bf83bd11c6bcc15f6d5e86a61f1ed9365fbaa0f4cf7"
    "gf256.cpp|fb9acc85b995714b0484fdfacf73427b3d4e33eef7c395f8745e876d01edfd90"
    "gf256.h|993b5de4a060b4efee4a45ac9ffc298e97b38c5efb00ead79dc21d87e8d470e2")
foreach(_pin IN LISTS _pins)
    string(REPLACE "|" ";" _fields "${_pin}")
    list(GET _fields 0 _file)
    list(GET _fields 1 _expected)
    file(SHA256 "${_root}/${_file}" _actual)
    if(NOT _actual STREQUAL _expected)
        message(FATAL_ERROR "Qualified input changed: ${_file}")
    endif()
endforeach()

# Reconstruct only the exact already-selected table; no scoring or codec calls.
execute_process(COMMAND /usr/bin/python3 "${_root}/bench/Wh2ThueMorseTinyGfniCostR1.py"
    lookup "${_binary}/lookup.bin" RESULT_VARIABLE _lookup_result TIMEOUT 30)
if(NOT _lookup_result STREQUAL "0")
    message(FATAL_ERROR "Fixed lookup reconstruction failed")
endif()
file(SHA256 "${_binary}/lookup.bin" _sha)
if(NOT _sha STREQUAL "27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d")
    message(FATAL_ERROR "Selected lookup identity differs")
endif()
file(READ "${_binary}/lookup.bin" _hex HEX)
string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," _literals "${_hex}")
file(WRITE "${_binary}/Wh2K6SerializedLookup.inc" "${_literals}\n")

foreach(_off IN ITEMS BUILD_SHARED_LIBS BUILD_TESTS BUILD_CODEC_V2 MARCH_NATIVE
        WIREHAIR_BUILD_BOTH WIREHAIR_BUILD_TOOLS WIREHAIR_BUILD_BENCHMARKS
        WIREHAIR_ENABLE_SCHEDULED_TESTS WIREHAIR_ENABLE_PUBLIC_BORROWED_R1_NATIVE_CONFIG
        WIREHAIR_ENABLE_THUE_MORSE_NATIVE_R0 WIREHAIR_ENABLE_LIBFUZZER
        WIREHAIR_ENABLE_TRI_BMQ4_SCHUR_FALSIFIER)
    set(${_off} OFF CACHE BOOL "" FORCE)
endforeach()
set(WIREHAIR_STATIC_PIC ON CACHE BOOL "" FORCE)
set(WIREHAIR_STRICT_WARNINGS ON CACHE BOOL "" FORCE)
set(WH_LTO OFF CACHE STRING "" FORCE)
set(WH_PGO_MODE OFF CACHE STRING "" FORCE)
add_subdirectory("${_root}" "${_binary}/production" EXCLUDE_FROM_ALL)

add_library(wh2_k6_serialized STATIC "${_root}/bench/Wh2K6Serialized.cpp"
    "${_root}/bench/Wh2ThueMorseTinyGfniR0Codec.cpp"
    "${_root}/bench/Wh2ThueMorseTinyGfniR0.cpp")
target_include_directories(wh2_k6_serialized PUBLIC "${_root}/bench")
target_link_libraries(wh2_k6_serialized PUBLIC Threads::Threads)
add_executable(wh2_k6_serialized_test "${_root}/bench/Wh2K6SerializedTest.cpp")
add_executable(wh2_k6_serialized_c_test "${_root}/bench/Wh2K6SerializedCTest.c")
set_target_properties(wh2_k6_serialized_c_test PROPERTIES LINKER_LANGUAGE CXX)
foreach(_target IN ITEMS wh2_k6_serialized_test wh2_k6_serialized_c_test)
    # The candidate runtime resolves GF before the public archive is scanned.
    # Only one GF runtime is linked; the archive's gf256.o remains unselected.
    target_link_libraries(${_target} PRIVATE wh2_k6_serialized wirehair)
    target_link_options(${_target} PRIVATE -fno-lto "-Wl,-Map,${_binary}/${_target}.map")
    add_test(NAME ${_target} COMMAND ${_target})
    set_tests_properties(${_target} PROPERTIES TIMEOUT 60)
    if(WH2_K6_SANITIZERS)
        target_link_options(${_target} PRIVATE -fsanitize=address,undefined -no-pie)
        set_tests_properties(${_target} PROPERTIES ENVIRONMENT
            "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:detect_stack_use_after_return=1;UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1")
    endif()
endforeach()
foreach(_target IN ITEMS wh2_k6_serialized wirehair wh2_k6_serialized_test wh2_k6_serialized_c_test)
    target_include_directories(${_target} PRIVATE "${_root}" "${_binary}")
    set_target_properties(${_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_compile_options(${_target} PRIVATE -Wall -Wextra -Wpedantic -Werror -fno-lto)
    if(WH2_K6_SCALAR)
        # Select the existing portable backend, not an Android-platform claim.
        target_compile_definitions(${_target} PRIVATE ANDROID WH2_K6_EXPECT_SCALAR=1)
    endif()
    if(WH2_K6_SANITIZERS)
        target_compile_options(${_target} PRIVATE -O1 -g -fno-omit-frame-pointer
            -fsanitize=address,undefined -fno-pie)
    else()
        target_compile_options(${_target} PRIVATE -O3)
    endif()
endforeach()
