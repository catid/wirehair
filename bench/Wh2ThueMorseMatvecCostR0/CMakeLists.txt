cmake_minimum_required(VERSION 3.20)
project(Wh2ThueMorseMatvecCostR0 LANGUAGES CXX)

option(WH2_MATVEC_COST_SANITIZERS "Use ASan/UBSan; scientific worker is disabled" OFF)
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
set(_flags "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_config}} ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_STATIC_LINKER_FLAGS}")
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

# Authenticate the reused construction and production source-list definitions.
# The neutral subproject independently checks all eight original codec/GF pins.
set(_pins
    "CMakeLists.txt|52742b1c4253185688f85340cd2be13eb82df590c7e6cacef520768c0477468a"
    "bench/Wh2ThueMorseMatvecCodecR0/CMakeLists.txt|e78a7ff61a9b67adf7995e6a2deb6ac23a6baa634114b54fb671b478010d7761"
    "bench/Wh2ThueMorseMatvecCodecR0Test.cpp|d4141849843b701b1faf6c1e4d9319043205eafcb80667a6752ad8ed6fe89c48")
foreach(_pin IN LISTS _pins)
    string(REPLACE "|" ";" _fields "${_pin}")
    list(GET _fields 0 _relative)
    list(GET _fields 1 _expected)
    execute_process(COMMAND /usr/bin/test -f "${_root}/${_relative}"
        RESULT_VARIABLE _regular TIMEOUT 1)
    if(IS_SYMLINK "${_root}/${_relative}" OR NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Frozen input must be regular: ${_relative}")
    endif()
    file(SHA256 "${_root}/${_relative}" _actual)
    if(NOT _actual STREQUAL _expected)
        message(FATAL_ERROR "Frozen input changed: ${_relative}")
    endif()
endforeach()

# Record the whole native repository dependency closure, plus build/controller
# sources, before generating or compiling anything. No old cost worker is used.
set(_native_inputs
    CMakeLists.txt cmake/wirehair.pc.in cmake/wirehairConfig.cmake.in
    wirehair.cpp include/wirehair/wirehair.h include/wirehair/wirehair.hpp
    gf256.cpp gf256.h WirehairCodec.cpp WirehairCodec.h WirehairCompiler.h
    WirehairEnvironment.h WirehairHeavy.h WirehairTools.cpp WirehairTools.h
    WirehairDenseFixups.inc WirehairPeelFixups.inc
    codec/WirehairV2Codec.cpp codec/WirehairV2Codec.h
    codec/WirehairV2Peel.cpp codec/WirehairV2Peel.h
    codec/WirehairV2Plan.cpp codec/WirehairV2Plan.h
    codec/WirehairV2Policy.cpp codec/WirehairV2Policy.h
    codec/WirehairV2Precode.cpp codec/WirehairV2Precode.h
    codec/WirehairV2PrecodeDecode.cpp codec/WirehairV2PrecodeDecode.h
    codec/WirehairV2PrecodeEncode.cpp codec/WirehairV2PrecodeEncode.h
    codec/WirehairV2Profile.cpp codec/WirehairV2Seeds.cpp codec/WirehairV2Seeds.h
    codec/WirehairV2Solve.cpp codec/WirehairV2Solve.h
    bench/Wh2FrozenTrace.cpp bench/Wh2FrozenTrace.h
    bench/Wh2PublicBorrowedTargetIdentity.cpp bench/Wh2PublicBorrowedTargetIdentity.h
    bench/Wh2RdpruTargetIdentityV2.cpp bench/Wh2RdpruTargetIdentityV2.h
    bench/Wh2ThueMorseNativeCodec.cpp bench/Wh2ThueMorseNativeCodec.h
    bench/Wh2ThueMorseNativeCodecTest.cpp
    bench/Wh2ThueMorseMatvecGfniR0.cpp bench/Wh2ThueMorseMatvecGfniR0.h
    bench/Wh2ThueMorseMatvecCodecR0/CMakeLists.txt bench/Wh2ThueMorseMatvecCodecR0Test.cpp
    bench/Wh2ThueMorseMatvecCostR0.cpp bench/Wh2ThueMorseMatvecCostR0Bridge.cpp
    bench/Wh2ThueMorseMatvecCostR0Bridge.h bench/Wh2ThueMorseMatvecCostR0BridgeTest.cpp
    bench/Wh2ThueMorseMatvecCostR0/CMakeLists.txt
    bench/Wh2ThueMorseMatvecCostR0.py bench/test_Wh2ThueMorseMatvecCostR0.py)
list(SORT _native_inputs)
set(_manifest "{\"schema\":\"wirehair.wh2.thue-matvec-cost-r0.native-inputs.v1\",\"sources\":{")
set(_separator "")
set(_input_hashes "")
foreach(_relative IN LISTS _native_inputs)
    set(_path "${_root}/${_relative}")
    execute_process(COMMAND /usr/bin/test -f "${_path}"
        RESULT_VARIABLE _regular TIMEOUT 1)
    if(IS_SYMLINK "${_path}" OR NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Native inputs must be regular source files: ${_relative}")
    endif()
    file(SHA256 "${_path}" _digest)
    list(APPEND _input_hashes "${_digest}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_path}")
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
        message(FATAL_ERROR "Native input manifest must be regular")
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
    if(EXISTS "${_binary}/build.ninja" OR EXISTS "${_binary}/Makefile")
        message(FATAL_ERROR "Cannot adopt an unrecorded existing build")
    endif()
    file(WRITE "${_manifest_path}" "${_manifest}")
endif()

# Reuse the qualified literal Apply substitution and its generated-source seal.
set(WH2_MATVEC_CODEC_SANITIZERS ${WH2_MATVEC_COST_SANITIZERS} CACHE BOOL "" FORCE)
add_subdirectory("${_root}/bench/Wh2ThueMorseMatvecCodecR0" "${_binary}/codec-neutral")
set(_generated "${_binary}/codec-neutral/Wh2ThueMorseMatvecCodecR0.generated.cpp")
file(SHA256 "${_generated}" _generated_sha)
if(NOT _generated_sha STREQUAL "2bc9b40321a2a95af982d66c7e9b3623f57674ba80ca3faf3494d9a54e9f84a8")
    message(FATAL_ERROR "Qualified generated codec differs")
endif()

# BUILD_CODEC_V2 disables internal experiments, not the public WH2 codec. The
# unchanged archive supplies both real public APIs. None of its old workers run.
foreach(_off IN ITEMS BUILD_SHARED_LIBS BUILD_TESTS BUILD_CODEC_V2 MARCH_NATIVE
        WIREHAIR_BUILD_BOTH WIREHAIR_BUILD_TOOLS WIREHAIR_BUILD_BENCHMARKS
        WIREHAIR_ENABLE_SCHEDULED_TESTS WIREHAIR_ENABLE_PUBLIC_BORROWED_R1_NATIVE_CONFIG
        WIREHAIR_ENABLE_THUE_MORSE_NATIVE_R0 WIREHAIR_ENABLE_LIBFUZZER
        WIREHAIR_ENABLE_TRI_BMQ4_SCHUR_FALSIFIER)
    set(${_off} OFF CACHE BOOL "Disabled in isolated matvec cost build" FORCE)
endforeach()
set(WIREHAIR_STATIC_PIC ON CACHE BOOL "" FORCE)
set(WIREHAIR_STRICT_WARNINGS ON CACHE BOOL "" FORCE)
set(WH_LTO OFF CACHE STRING "" FORCE)
set(WH_PGO_MODE OFF CACHE STRING "" FORCE)
add_subdirectory("${_root}" "${_binary}/production" EXCLUDE_FROM_ALL)
# The unchanged project defines these five offline generators even with tools
# disabled. They are not built here and must not enter the actual-TU receipt.
# This target-scoped property requires CMake 3.20; no compile/link flag changes.
foreach(_generator IN ITEMS gen_small_dseeds gen_peel_seeds gen_most_dseeds gen_dcounts gen_tables)
    if(NOT TARGET ${_generator})
        message(FATAL_ERROR "Expected unbuilt offline generator is absent")
    endif()
    get_target_property(_excluded ${_generator} EXCLUDE_FROM_ALL)
    if(NOT _excluded)
        message(FATAL_ERROR "Offline generator unexpectedly participates in the build")
    endif()
    set_target_properties(${_generator} PROPERTIES EXPORT_COMPILE_COMMANDS OFF)
endforeach()

# Reconstruction is confined to the already selected lookup. It is not a codec
# run, a parameter selection, or a correctness/performance campaign.
execute_process(COMMAND /usr/bin/python3 "${_root}/bench/Wh2ThueMorseMatvecCostR0.py"
    lookup "${_binary}/lookup.bin" RESULT_VARIABLE _lookup_result TIMEOUT 30)
if(NOT _lookup_result STREQUAL "0")
    message(FATAL_ERROR "Fixed lookup reconstruction failed")
endif()
execute_process(COMMAND /usr/bin/test -f "${_binary}/lookup.bin"
    RESULT_VARIABLE _lookup_regular TIMEOUT 1)
if(IS_SYMLINK "${_binary}/lookup.bin" OR NOT _lookup_regular STREQUAL "0")
    message(FATAL_ERROR "Lookup must be a regular non-symlink file")
endif()
file(SIZE "${_binary}/lookup.bin" _lookup_size)
file(SHA256 "${_binary}/lookup.bin" _lookup_sha)
if(NOT _lookup_size EQUAL 39936 OR
   NOT _lookup_sha STREQUAL "27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d")
    message(FATAL_ERROR "Selected lookup bytes differ")
endif()

set(_bridge_objects "")
set(_bridges "")
foreach(_arm RANGE 0 3)
    set(_target "matvec_cost_bridge_${_arm}")
    add_library(${_target} OBJECT "${_root}/bench/Wh2ThueMorseMatvecCostR0Bridge.cpp")
    target_compile_definitions(${_target} PRIVATE WH2_MATVEC_COST_ARM=${_arm})
    target_link_libraries(${_target} PRIVATE wirehair)
    list(APPEND _bridges ${_target})
    list(APPEND _bridge_objects $<TARGET_OBJECTS:${_target}>)
endforeach()
set(_codec_objects $<TARGET_OBJECTS:matvec_codec_baseline> $<TARGET_OBJECTS:matvec_codec_candidate>)
add_executable(wh2_thue_matvec_cost_r0
    "${_root}/bench/Wh2ThueMorseMatvecCostR0.cpp"
    "${_root}/bench/Wh2RdpruTargetIdentityV2.cpp"
    "${_root}/bench/Wh2PublicBorrowedTargetIdentity.cpp"
    "${_root}/bench/Wh2FrozenTrace.cpp" ${_codec_objects} ${_bridge_objects})
target_compile_definitions(wh2_thue_matvec_cost_r0 PRIVATE
    WH2_MATVEC_COST_LOOKUP_PATH="${_binary}/lookup.bin"
    WH2_MATVEC_COST_SANITIZERS=$<BOOL:${WH2_MATVEC_COST_SANITIZERS}>)
target_link_options(wh2_thue_matvec_cost_r0 PRIVATE "-Wl,-Map,${_binary}/matvec-cost.map")
add_executable(wh2_matvec_cost_bridge_test
    "${_root}/bench/Wh2ThueMorseMatvecCostR0BridgeTest.cpp" ${_codec_objects} ${_bridge_objects})
target_link_options(wh2_matvec_cost_bridge_test PRIVATE "-Wl,-Map,${_binary}/matvec-cost-bridge-test.map")
set(_executables wh2_thue_matvec_cost_r0 wh2_matvec_cost_bridge_test)
foreach(_target IN ITEMS wirehair ${_bridges} ${_executables})
    target_include_directories(${_target} PRIVATE "${_root}" "${_root}/bench" "${_root}/include")
    set_target_properties(${_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_compile_options(${_target} PRIVATE -Wall -Wextra -Wpedantic -Werror -fno-lto)
    if(WH2_MATVEC_COST_SANITIZERS)
        target_compile_options(${_target} PRIVATE -O1 -g -fno-omit-frame-pointer
            -fsanitize=address,undefined -fno-pie)
    else()
        target_compile_options(${_target} PRIVATE -O3)
    endif()
endforeach()
foreach(_target IN LISTS _executables)
    # The .72 archive is searched first and supplies all GF symbols. The later
    # production archive's gf256.cpp.o must remain unextracted (audit link map).
    target_link_libraries(${_target} PRIVATE matvec_codec_gf wirehair Threads::Threads)
    target_link_options(${_target} PRIVATE -fno-lto)
    if(WH2_MATVEC_COST_SANITIZERS)
        target_link_options(${_target} PRIVATE -fsanitize=address,undefined -no-pie)
    endif()
endforeach()
add_test(NAME matvec_cost_neutral COMMAND wh2_thue_matvec_cost_r0 --neutral)
add_test(NAME matvec_cost_bridge_neutral COMMAND wh2_matvec_cost_bridge_test --require-gfni)
set_tests_properties(matvec_cost_neutral matvec_cost_bridge_neutral PROPERTIES TIMEOUT 60)
if(WH2_MATVEC_COST_SANITIZERS)
    set_tests_properties(matvec_cost_neutral matvec_cost_bridge_neutral PROPERTIES ENVIRONMENT
        "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:detect_stack_use_after_return=1;UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1")
else()
    find_program(_cost_qemu qemu-x86_64)
    if(_cost_qemu)
        add_test(NAME matvec_cost_bridge_haswell COMMAND "${_cost_qemu}" -cpu Haswell
            $<TARGET_FILE:wh2_matvec_cost_bridge_test> --expect-scalar)
        set_tests_properties(matvec_cost_bridge_haswell PROPERTIES TIMEOUT 60)
    endif()
endif()

# A configure-time edit must not produce a manifest for one set of bytes and
# compilation rules for another. The controller additionally checks depfiles.
set(_index 0)
foreach(_relative IN LISTS _native_inputs)
    list(GET _input_hashes ${_index} _before)
    execute_process(COMMAND /usr/bin/test -f "${_root}/${_relative}"
        RESULT_VARIABLE _regular TIMEOUT 1)
    if(IS_SYMLINK "${_root}/${_relative}" OR NOT _regular STREQUAL "0")
        message(FATAL_ERROR "Native input changed type during configuration")
    endif()
    file(SHA256 "${_root}/${_relative}" _after)
    if(NOT _before STREQUAL _after)
        message(FATAL_ERROR "Native input changed during configuration: ${_relative}")
    endif()
    math(EXPR _index "${_index} + 1")
endforeach()
