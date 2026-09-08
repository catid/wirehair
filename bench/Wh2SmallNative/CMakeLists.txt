cmake_minimum_required(VERSION 3.16)
project(Wh2SmallNative LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
get_filename_component(WH2_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(WH2_SMALL_LIBRARY "" CACHE FILEPATH "Exact existing qualified Wirehair static library")
option(WH2_SMALL_PORTABLE "Match an ANDROID-selected portable GF256 library layout (not an Android platform test)" OFF)
if(NOT EXISTS "${WH2_SMALL_LIBRARY}")
    message(FATAL_ERROR "WH2_SMALL_LIBRARY must name an existing qualified library")
endif()
find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter)
add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/Wh2K3NativeData.inc"
    COMMAND "${Python3_EXECUTABLE}" -I -B -S "${WH2_ROOT}/bench/Wh2K3NativeDataR0.py"
            --header "${CMAKE_CURRENT_BINARY_DIR}/Wh2K3NativeData.inc"
    # The generator preserves matching bytes. Refresh only the build output's
    # timestamp so an unchanged reconstruction does not run on every build.
    COMMAND "${CMAKE_COMMAND}" -E touch "${CMAKE_CURRENT_BINARY_DIR}/Wh2K3NativeData.inc"
    DEPENDS "${WH2_ROOT}/bench/Wh2K3NativeDataR0.py"
            "${WH2_ROOT}/bench/Wh2NoncommutingRadixRunR0.py"
    VERBATIM)
add_executable(wh2_small_native_test
    "${WH2_ROOT}/bench/Wh2SmallNativeTest.cpp"
    "${WH2_ROOT}/bench/Wh2FrozenTrace.cpp"
    "${CMAKE_CURRENT_BINARY_DIR}/Wh2K3NativeData.inc")
target_include_directories(wh2_small_native_test PRIVATE "${WH2_ROOT}" "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(wh2_small_native_test PRIVATE "${WH2_SMALL_LIBRARY}")
if(WH2_SMALL_PORTABLE)
    target_compile_definitions(wh2_small_native_test PRIVATE ANDROID WH2_SMALL_EXPECT_PORTABLE=1)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(wh2_small_native_test PRIVATE -Wall -Wextra -Wpedantic -Werror -fno-strict-aliasing)
endif()
enable_testing()
add_test(NAME small_neutral COMMAND wh2_small_native_test --neutral)
add_test(NAME small_corpus COMMAND wh2_small_native_test --corpus)
set_tests_properties(small_neutral small_corpus PROPERTIES TIMEOUT 120)
