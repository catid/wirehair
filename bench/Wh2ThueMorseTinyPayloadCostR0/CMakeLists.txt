cmake_minimum_required(VERSION 3.16)
project(Wh2ThueMorseTinyPayloadCostR0 LANGUAGES CXX)

option(WH2_TINY_COST_SANITIZERS "Use ASan/UBSan; scientific worker is then disabled" OFF)
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
   _flags MATCHES "WH_COUNT|WIREHAIR_TESTING|WIREHAIR_V2_ENABLE_TEST_HOOKS|-flto|-fprofile|-march=|-mavx|-mgfni")
    message(FATAL_ERROR "Counters, hooks, LTO, PGO and injected ISA flags are forbidden")
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(Threads REQUIRED)
enable_testing()

# Immutable prebuild native inputs: a changed source requires a fresh build,
# even if its mtime was preserved. Never adopt an existing unrecorded build.
set(_native_inputs gf256.cpp gf256.h WirehairEnvironment.h
    bench/Wh2FrozenTrace.cpp bench/Wh2FrozenTrace.h
    bench/Wh2PublicBorrowedTargetIdentity.cpp bench/Wh2PublicBorrowedTargetIdentity.h
    bench/Wh2RdpruTargetIdentityV2.cpp bench/Wh2RdpruTargetIdentityV2.h
    bench/Wh2ThueMorseNativeCodec.cpp bench/Wh2ThueMorseNativeCodec.h
    bench/Wh2ThueMorseNativeCodecTest.cpp bench/Wh2ThueMorseTinyPayloadR0.cpp
    bench/Wh2ThueMorseTinyPayloadR0.h bench/Wh2ThueMorseTinyPayloadR0CodecTest.cpp
    bench/Wh2ThueMorseTinyPayloadR0Test.cpp bench/Wh2ThueMorseTinyPayloadCostR0.cpp
    bench/Wh2ThueMorseTinyPayloadCostR0Bridge.cpp bench/Wh2ThueMorseTinyPayloadCostR0Bridge.h)
list(SORT _native_inputs)
set(_manifest "{\"schema\":\"wirehair.wh2.thue-tiny-payload-cost-r0.native-inputs.v1\",\"sources\":{")
set(_separator "")
foreach(_relative IN LISTS _native_inputs)
    execute_process(COMMAND /usr/bin/test -f "${_root}/${_relative}"
        RESULT_VARIABLE _regular TIMEOUT 1)
    if(IS_SYMLINK "${_root}/${_relative}" OR NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Native inputs must be regular source files")
    endif()
    file(SHA256 "${_root}/${_relative}" _digest)
    string(APPEND _manifest "${_separator}\"${_relative}\":\"${_digest}\"")
    set(_separator ",")
endforeach()
string(APPEND _manifest "}}\n")
set(_manifest_path "${_binary}/native-inputs.json")
if(IS_SYMLINK "${_manifest_path}" OR IS_DIRECTORY "${_manifest_path}")
    message(FATAL_ERROR "Invalid native input manifest")
endif()
if(EXISTS "${_manifest_path}")
    execute_process(COMMAND /usr/bin/test -f "${_manifest_path}"
        RESULT_VARIABLE _regular TIMEOUT 1)
    if(NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Native input manifest must be a regular file")
    endif()
    file(SIZE "${_manifest_path}" _manifest_size)
    if(_manifest_size GREATER 65536)
        message(FATAL_ERROR "Oversized native input manifest")
    endif()
    file(READ "${_manifest_path}" _old_manifest)
    if(NOT _old_manifest STREQUAL _manifest)
        message(FATAL_ERROR "Native inputs changed: use a fresh build directory")
    endif()
else()
    if(EXISTS "${_binary}/build.ninja")
        message(FATAL_ERROR "Cannot adopt an unrecorded existing build")
    endif()
    file(WRITE "${_manifest_path}" "${_manifest}")
endif()

# Reconstruct only the already selected table. This never executes a codec,
# selects a parameter, samples loss, or evaluates recovery/performance.
set(_controller "${_root}/bench/Wh2ThueMorseTinyPayloadCostR0.py")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_controller}")
execute_process(COMMAND /usr/bin/python3 "${_controller}" lookup "${_binary}/lookup.bin"
    RESULT_VARIABLE _lookup_result TIMEOUT 30)
if(NOT _lookup_result STREQUAL "0")
    message(FATAL_ERROR "Fixed lookup reconstruction failed")
endif()

add_library(tiny_cost_gf STATIC "${_root}/gf256.cpp")
target_link_libraries(tiny_cost_gf PUBLIC Threads::Threads)
add_library(tiny_cost_baseline OBJECT "${_root}/bench/Wh2ThueMorseNativeCodec.cpp")
add_library(tiny_cost_candidate OBJECT "${_root}/bench/Wh2ThueMorseTinyPayloadR0.cpp")
add_library(tiny_cost_baseline_bridge OBJECT "${_root}/bench/Wh2ThueMorseTinyPayloadCostR0Bridge.cpp")
target_compile_definitions(tiny_cost_baseline_bridge PRIVATE WH2_TINY_COST_CANDIDATE=0)
add_library(tiny_cost_candidate_bridge OBJECT "${_root}/bench/Wh2ThueMorseTinyPayloadCostR0Bridge.cpp")
target_compile_definitions(tiny_cost_candidate_bridge PRIVATE WH2_TINY_COST_CANDIDATE=1)
add_executable(wh2_thue_tiny_payload_cost_r0
    "${_root}/bench/Wh2ThueMorseTinyPayloadCostR0.cpp"
    "${_root}/bench/Wh2RdpruTargetIdentityV2.cpp"
    "${_root}/bench/Wh2PublicBorrowedTargetIdentity.cpp"
    "${_root}/bench/Wh2FrozenTrace.cpp"
    $<TARGET_OBJECTS:tiny_cost_baseline> $<TARGET_OBJECTS:tiny_cost_candidate>
    $<TARGET_OBJECTS:tiny_cost_baseline_bridge> $<TARGET_OBJECTS:tiny_cost_candidate_bridge>)
target_compile_definitions(wh2_thue_tiny_payload_cost_r0 PRIVATE
    WH2_TINY_COST_LOOKUP_PATH="${_binary}/lookup.bin"
    WH2_TINY_COST_SANITIZERS=$<BOOL:${WH2_TINY_COST_SANITIZERS}>)
target_link_options(wh2_thue_tiny_payload_cost_r0 PRIVATE "-Wl,-Map,${_binary}/tiny-cost.map")
add_executable(tiny_cost_baseline_test "${_root}/bench/Wh2ThueMorseNativeCodecTest.cpp"
    $<TARGET_OBJECTS:tiny_cost_baseline>)
add_executable(tiny_cost_candidate_test "${_root}/bench/Wh2ThueMorseTinyPayloadR0CodecTest.cpp"
    $<TARGET_OBJECTS:tiny_cost_candidate>)
add_executable(tiny_cost_kernel_test "${_root}/bench/Wh2ThueMorseTinyPayloadR0Test.cpp")
set(_executables wh2_thue_tiny_payload_cost_r0 tiny_cost_baseline_test tiny_cost_candidate_test tiny_cost_kernel_test)
foreach(_target IN ITEMS tiny_cost_gf tiny_cost_baseline tiny_cost_candidate
        tiny_cost_baseline_bridge tiny_cost_candidate_bridge ${_executables})
    target_include_directories(${_target} PRIVATE "${_root}" "${_root}/bench" "${_root}/include")
    set_target_properties(${_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_compile_options(${_target} PRIVATE -Wall -Wextra -Wpedantic -Werror -fno-lto)
    if(WH2_TINY_COST_SANITIZERS)
        target_compile_options(${_target} PRIVATE -O1 -g -fno-omit-frame-pointer
            -fsanitize=address,undefined -fno-pie)
    else()
        target_compile_options(${_target} PRIVATE -O3)
    endif()
endforeach()
foreach(_target IN LISTS _executables)
    target_link_libraries(${_target} PRIVATE tiny_cost_gf)
    target_link_options(${_target} PRIVATE -fno-lto)
    if(WH2_TINY_COST_SANITIZERS)
        target_link_options(${_target} PRIVATE -fsanitize=address,undefined -no-pie)
    endif()
endforeach()
add_test(NAME tiny_cost_neutral COMMAND wh2_thue_tiny_payload_cost_r0 --neutral)
foreach(_target IN ITEMS tiny_cost_baseline_test tiny_cost_candidate_test tiny_cost_kernel_test)
    add_test(NAME ${_target} COMMAND ${_target})
endforeach()
set_tests_properties(tiny_cost_neutral tiny_cost_baseline_test tiny_cost_candidate_test tiny_cost_kernel_test
    PROPERTIES TIMEOUT 60)
