cmake_minimum_required(VERSION 3.16)
project(Wh2SmallNative LANGUAGES C CXX)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
get_filename_component(WH2_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(WH2_SMALL_LIBRARY "" CACHE FILEPATH "Exact existing qualified Wirehair static library")
option(WH2_SMALL_PORTABLE "Match an ANDROID-selected portable GF256 library layout (not an Android platform test)" OFF)
if(NOT EXISTS "${WH2_SMALL_LIBRARY}")
    message(FATAL_ERROR "WH2_SMALL_LIBRARY must name an existing qualified library")
endif()
set(WH2_SMALL_TEST_DIMENSION "3" CACHE STRING "Sealed benchmark dimension: 2, 3, 4, 5 or 8")
set_property(CACHE WH2_SMALL_TEST_DIMENSION PROPERTY STRINGS 2 3 4 5 8)
if(NOT WH2_SMALL_TEST_DIMENSION MATCHES "^(2|3|4|5|8)$")
    message(FATAL_ERROR "Only sealed benchmark dimensions 2, 3, 4, 5 and 8 are supported")
endif()
set(WH2_SMALL_DATA "Wh2K${WH2_SMALL_TEST_DIMENSION}NativeData")
find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter)
add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${WH2_SMALL_DATA}.inc"
    COMMAND "${Python3_EXECUTABLE}" -I -B -S "${WH2_ROOT}/bench/${WH2_SMALL_DATA}R0.py"
            --header "${CMAKE_CURRENT_BINARY_DIR}/${WH2_SMALL_DATA}.inc"
    # The generator preserves matching bytes. Refresh only the build output's
    # timestamp so an unchanged reconstruction does not run on every build.
    COMMAND "${CMAKE_COMMAND}" -E touch "${CMAKE_CURRENT_BINARY_DIR}/${WH2_SMALL_DATA}.inc"
    DEPENDS "${WH2_ROOT}/bench/${WH2_SMALL_DATA}R0.py"
            "${WH2_ROOT}/bench/Wh2NoncommutingRadixRunR0.py"
            "${WH2_ROOT}/bench/Wh2K3NativeDataR0.py"
    VERBATIM)
add_custom_target(wh2_small_data DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${WH2_SMALL_DATA}.inc")
add_executable(wh2_small_native_test
    "${WH2_ROOT}/bench/Wh2SmallNativeTest.cpp"
    "${WH2_ROOT}/bench/Wh2FrozenTrace.cpp")
add_dependencies(wh2_small_native_test wh2_small_data)
target_include_directories(wh2_small_native_test PRIVATE "${WH2_ROOT}" "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(wh2_small_native_test PRIVATE "${WH2_SMALL_LIBRARY}")
if(WH2_SMALL_PORTABLE)
    target_compile_definitions(wh2_small_native_test PRIVATE ANDROID WH2_SMALL_EXPECT_PORTABLE=1)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(wh2_small_native_test PRIVATE -Wall -Wextra -Wpedantic -Werror -fno-strict-aliasing)
endif()
target_compile_definitions(wh2_small_native_test PRIVATE WH2_SMALL_CODEC_K=${WH2_SMALL_TEST_DIMENSION})
enable_testing()
add_test(NAME small_neutral COMMAND wh2_small_native_test --neutral)
add_test(NAME small_corpus COMMAND wh2_small_native_test --corpus)
set_tests_properties(small_neutral small_corpus PROPERTIES TIMEOUT 120)

# Keep the complete serialized C boundary in a separate, non-LTO translation
# unit. The installed library and its qualified K6 implementation are inputs,
# not edited/rebuilt/replaced by this standalone project.
add_library(wh2_small_serialized STATIC "${WH2_ROOT}/bench/Wh2SmallSerialized.cpp")
add_dependencies(wh2_small_serialized wh2_small_data)
set_target_properties(wh2_small_serialized PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(wh2_small_serialized PUBLIC "${WH2_SMALL_LIBRARY}")
add_executable(wh2_small_serialized_test "${WH2_ROOT}/bench/Wh2SmallSerializedTest.cpp")
add_executable(wh2_small_serialized_k6_test "${WH2_ROOT}/bench/Wh2SmallSerializedTest.cpp")
target_compile_definitions(wh2_small_serialized_k6_test PRIVATE WH2_SMALL_TEST_K=6)
add_executable(wh2_small_serialized_c_test "${WH2_ROOT}/bench/Wh2SmallSerializedCTest.c")
set_target_properties(wh2_small_serialized_c_test PROPERTIES LINKER_LANGUAGE CXX)
add_executable(wh2_small_serialized_corpus_test
    "${WH2_ROOT}/bench/Wh2SmallNativeTest.cpp" "${WH2_ROOT}/bench/Wh2FrozenTrace.cpp")
add_dependencies(wh2_small_serialized_corpus_test wh2_small_data)
target_compile_definitions(wh2_small_serialized_corpus_test PRIVATE WH2_SMALL_TEST_SERIALIZED=1)
foreach(target IN ITEMS wh2_small_serialized wh2_small_serialized_test wh2_small_serialized_k6_test
        wh2_small_serialized_c_test wh2_small_serialized_corpus_test)
    target_compile_definitions(${target} PRIVATE WH2_SMALL_CODEC_K=${WH2_SMALL_TEST_DIMENSION})
    target_include_directories(${target} PRIVATE "${WH2_ROOT}" "${WH2_ROOT}/include" "${CMAKE_CURRENT_BINARY_DIR}")
    if(WH2_SMALL_PORTABLE)
        target_compile_definitions(${target} PRIVATE ANDROID WH2_SMALL_EXPECT_PORTABLE=1)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror -fno-strict-aliasing -fno-lto)
        target_link_options(${target} PRIVATE -fno-lto)
    endif()
    if(NOT target STREQUAL "wh2_small_serialized")
        target_link_libraries(${target} PRIVATE wh2_small_serialized)
    endif()
endforeach()
add_test(NAME small_serialized COMMAND wh2_small_serialized_test)
add_test(NAME small_serialized_k6 COMMAND wh2_small_serialized_k6_test)
add_test(NAME small_serialized_c COMMAND wh2_small_serialized_c_test)
add_test(NAME small_serialized_neutral COMMAND wh2_small_serialized_corpus_test --neutral)
add_test(NAME small_serialized_corpus COMMAND wh2_small_serialized_corpus_test --corpus)
set_tests_properties(small_serialized small_serialized_k6 small_serialized_c
    small_serialized_neutral small_serialized_corpus PROPERTIES TIMEOUT 120)
